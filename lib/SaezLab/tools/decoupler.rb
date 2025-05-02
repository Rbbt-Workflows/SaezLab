require 'rbbt-util'
require 'rbbt/util/python'
require 'rbbt/util/python/util'

module SaezLab

  def self._pycall_decoupler(&block)
    RbbtPython.run "decoupler", as: :decoupler_python, &block
  end
  
  def self.decoupler_methods
    %w(ulm mlm wsum combined)
  end

  def self.run_decoupler(method, mat, net, options = {})
    mat = mat.tsv if Path === mat
    net = net.tsv if Path === net

    mat = RbbtPython.tsv2df(mat) if TSV === mat
    net = RbbtPython.tsv2df(net) if TSV === net

    method = method.to_s
    method = 'decouple' if method.include?('combined')
    method = "run_" + method.to_s unless method =~ /^run_|decouple/ 

    acts = norm_acts = pvals = nil
    _pycall_decoupler do
      Log.low "Running RbbtPython code"
      mat = mat.astype('float')

      net["weight"] = net["weight"].astype('float') if RbbtPython.to_a(net.columns.values).include?('weight')

      Log.low "Calling decouple"
      case method
      when "run_mlm"
        res_tuple = decoupler_python.run_mlm(mat, net, **options)
      when "run_ulm"
        res_tuple = decoupler_python.run_ulm(mat, net, **options)
      when "decouple"
        res_tuple = decoupler_python.decouple(mat, net, **options)
      else
        res_tuple = decoupler_python.send(method.to_sym, mat, net, **options)
      end
      Log.low "Done decouple"

      if method == 'decouple'
        acts = nil
        res_tuple.each do |name, d|
          new = df2tsv(d, :cast => :to_f) 
          new.fields = new.fields.collect{|f| "#{f} (#{name})" }
          if acts.nil?
            acts = new
          else
            acts.attach new
          end
        end
        norm_acts = pvals = nil
      else

        len = PyCall.len(res_tuple)
        acts = df2tsv(res_tuple[0], :cast => :to_f) 
        if len == 2
          pvals = df2tsv(res_tuple[1], :cast => :to_f)
        else
          norm_acts = df2tsv(res_tuple[1], :cast => :to_f)
          pvals = df2tsv(res_tuple[2], :cast => :to_f)
        end
      end
    end

    [acts, norm_acts, pvals]
  end

  def self.run_decoupler_script(method, mat, net, options = {})
    mat = mat.tsv if Path === mat
    net = net.tsv if Path === net

    method = method.to_s
    method = 'decouple' if method.include?('combined')
    method = "run_" + method.to_s unless method =~ /^run_|decouple/ 

    TmpFile.with_dir do |dir|
      Path.setup(dir)
      acts = norm_acts = pvals = nil
      RbbtPython.script <<-EOF, mat: mat, net: net, method: method, options: options, dir: dir
import decoupler

mat = mat.astype('float')

#res_tuple = decoupler.run_ulm(mat, net, **options)
res_tuple = getattr(decoupler, method)(mat, net, **options)

i = 0
for t in res_tuple:
  i = i + 1
  filename = dir + "/" + "file" + str(i) + ".tsv"
  rbbt.save_tsv(filename, t)
      EOF
      acts, norm_acts, pvals = dir.glob("*.tsv").collect{|f| TSV.open(f, type: :list, cast: :to_f) }
      pvals, norm_acts = norm_acts, nil if pvals.nil?
      [acts, norm_acts, pvals]
    end
  end
end
