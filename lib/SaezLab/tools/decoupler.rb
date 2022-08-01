require 'rbbt-util'
require 'rbbt/util/python'
require 'rbbt/util/python/util'

module SaezLab

  def self._pycall_decoupler(&block)
    RbbtPython.init_rbbt
    RbbtPython.run "decoupler" do 
      RbbtPython.module_eval(&block)
    end
  end
  
  def self.decoupler_methods
    @@decoupler_methods ||= Persist.persist("Decouple methods", :array) do
      _pycall_decoupler do
        decoupler.show_methods()["Function"].to_list()
      end
    end
  end

  def self.decoupler(method, mat, net, options = {})

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
      res_tuple = decoupler.send(method.to_sym, mat, net, **options)
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
end
