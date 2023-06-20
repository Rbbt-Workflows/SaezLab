module SaezLab

  dep :evaluate_knocktf, :negative_evidence_proportion => :sweep,
    :sign_support_balance => :sweep, :strict_negative => :sweep, :support_evidence_max => :sweep do |jobname,options|

    sign_support_balance_sweep = [0.0, 0.2, 0.3, 0.5, 0.7]
    support_evidence_max_sweep = [1,2,3]
    sign_support_evidence_max_sweep = [1,2,3]
    strict_negative_sweep = [0.2,0.3,0.4] 
    min_n_sweep = [3,4,5]
    negative_evidence_proportion_sweep = [1]
    use_tf_role_sweep = [true, false]

    jobs = []

    sign_support_balance_sweep.each do |sign_support_balance|
      strict_negative_sweep.each do |strict_negative|
        support_evidence_max_sweep.each do |support_evidence_max|
          sign_support_evidence_max_sweep.each do |sign_support_evidence_max|
            negative_evidence_proportion_sweep.each do |negative_evidence_proportion|
              min_n_sweep.each do |min_n|
                use_tf_role_sweep.each do |use_tf_role|

                  job_options = options.merge(
                    :sign_support_balance => sign_support_balance,
                    :strict_negative => strict_negative,
                    :support_evidence_max => support_evidence_max,
                    :sign_support_evidence_max => sign_support_evidence_max,
                    :negative_evidence_proportion => negative_evidence_proportion,
                    :use_tf_role => use_tf_role,
                    :min_n => min_n
                  )

                  jobs << {:inputs => job_options}
                end
              end
            end
          end
        end
      end
    end

    jobs
  end
  task :sweep => :tsv do

    fields = [:strict_negative, :negative_evidence_proportion, :use_tf_role, :support_evidence_max , :sign_support_evidence_max, :sign_support_balance, :min_n]
    dependencies.inject(nil) do |acc,dep|
      expression = fields.zip(dep.recursive_inputs.values_at *fields).collect{|f,v| [f,v] * "="} * ","
      tsv = dep.load.transpose("Experiment")
      tsv[expression] = tsv.delete "Value"
      
      acc = acc.nil? ? tsv : acc.merge!(tsv)
    end
  end


  dep :evaluate_knocktf_auroc, :negative_evidence_proportion => :sweep,
    :sign_support_balance => :sweep, :strict_negative => :sweep, :support_evidence_max => :sweep do |jobname,options|

    sign_support_balance_sweep = [0.0, 0.2, 0.3, 0.5, 0.7]
    support_evidence_max_sweep = [1,2,3]
    sign_support_evidence_max_sweep = [1,2,3]
    strict_negative_sweep = [0.2,0.3,0.4] 
    min_n_sweep = [3,4,5]
    negative_evidence_proportion_sweep = [1]
    use_tf_role_sweep = [true, false]

    jobs = []

    sign_support_balance_sweep.each do |sign_support_balance|
      strict_negative_sweep.each do |strict_negative|
        support_evidence_max_sweep.each do |support_evidence_max|
          sign_support_evidence_max_sweep.each do |sign_support_evidence_max|
            negative_evidence_proportion_sweep.each do |negative_evidence_proportion|
              min_n_sweep.each do |min_n|
                use_tf_role_sweep.each do |use_tf_role|

                  job_options = options.merge(
                    :sign_support_balance => sign_support_balance,
                    :strict_negative => strict_negative,
                    :support_evidence_max => support_evidence_max,
                    :sign_support_evidence_max => sign_support_evidence_max,
                    :negative_evidence_proportion => negative_evidence_proportion,
                    :use_tf_role => use_tf_role,
                    :min_n => min_n
                  )

                  jobs << {:inputs => job_options}
                end
              end
            end
          end
        end
      end
    end

    jobs
    end
  task :sweep_auroc => :tsv do

    fields = [:strict_negative, :negative_evidence_proportion, :use_tf_role, :support_evidence_max , :sign_support_evidence_max, :sign_support_balance, :min_n]
    tsv = TSV.setup({}, "Experimen~AUROC")
    dependencies.inject(tsv) do |acc,dep|
      expression = fields.zip(dep.recursive_inputs.values_at *fields).collect{|f,v| [f,v] * "="} * ","
      value = dep.load
      acc[expression] = value
    end
  end

  def self.ExTRI_databases
    Workflow.require_workflow "ExTRI"
    extri = ExTRI.job('pairs_final').load
    databases = extri.fields.collect do |f,v|
      m = f.match(/\[(.*)\] (.*)/)
      next unless m
      db, field = m.values_at 1, 2
      db
    end.flatten.uniq.compact
  end

  #{{{ --------------------- DATABASE SWEEPS ---------------------------
  
  #{{{ ++++ DATABASE SWEEPS +++

  dep :evaluate_knocktf, :databases => [], :jobname => "CollecTRI", :canfail => true
  dep :evaluate_knocktf, :databases => :placeholder, :canfail => true do |jobname,options|
    SaezLab.ExTRI_databases.collect do |db|
      {:inputs => options.merge(:databases => [db]), :jobname => db }
    end
  end
  dep :evaluate_knocktf, :databases => [], :dorothea => true, :jobname => "Dorothea"
  task :database_sweep => :tsv do
    dependencies.inject(nil) do |acc,dep|
      next acc if dep.error?
      if dep.recursive_inputs[:dorothea]
        database = "Dorothea"
      else
        database = dep.recursive_inputs[:databases].first
        if database.nil? || database.empty?
          database = "CollecTRI" 
        end
      end
      tsv = dep.load.transpose("Database")
      tsv[database] = tsv.delete "Value"
      
      acc = acc.nil? ? tsv : acc.merge!(tsv)
    end
  end

  dep :evaluate_knocktf, :databases => [], :jobname => "CollecTRI"
  dep :evaluate_knocktf, :databases => :placeholder, :canfail => true do |jobname,options|
    databases = SaezLab.ExTRI_databases
    databases.collect do |db|
      {:inputs => options.merge(:databases => databases - [db, "ExTRI", "CollecTRI", "Dorothea"]), :jobname => "Remove #{db}" }
    end
  end
  task :database_remove_sweep => :tsv do
    databases = SaezLab.ExTRI_databases

    dependencies.inject(nil) do |acc,dep|
      next acc if dep.error?
      database = (databases - dep.recursive_inputs[:databases]).first

      tsv = dep.load.transpose("Database")
      tsv[database] = tsv.delete "Value"
      
      acc = acc.nil? ? tsv : acc.merge!(tsv)
    end
  end

  #{{{ ++++ SWEEP BOXPLOT +++

  dep :database_sweep
  extension :png
  task :database_sweep_plot => :binary do
    R::PNG.ggplot self.tmp_path, step(:database_sweep).load, <<-EOR, 6, 4 
rbbt.require('colorspace')
data$Database = rownames(data)
names(data) <- make.names(names(data))
ggplot(data, ) + geom_point(aes(x=Matches,y=Accuracy..,color=Database)) +
   geom_text(aes(x=Matches,y=Accuracy..,label=Database)) + 
   theme_classic() +
   scale_fill_continuous_diverging()
    EOR
    nil
  end

  dep :database_remove_sweep
  extension :png
  task :database_remove_sweep_plot => :binary do
    R::PNG.ggplot self.tmp_path, step(:database_remove_sweep).load, <<-EOR, 6, 4
rbbt.require('colorspace')
data$Database = rownames(data)
names(data) <- make.names(names(data))
ggplot(data, ) + geom_point(aes(x=Matches,y=Accuracy..,color=Database)) +
   geom_text(aes(x=Matches,y=Accuracy..,label=Database)) + 
   theme_classic() +
   scale_fill_continuous_diverging()
    EOR
    nil
  end




  #{{{ ++++ R_BENCHMARK +++
   
  dep :regulome, :databases => [], :jobname => "CollecTRI"
  dep :regulome, :databases => :regulome, :compute => :produce do |jobname,options|
    SaezLab.ExTRI_databases.collect do |db|
      {:inputs => options.merge(:databases => [db]), :jobname => db }
    end
  end
  dep :dorothea, :jobname => "Dorothea"
  dep :decouple_R_benchmark, :network => :sweep, :canfail => true do |jobname,options,dependencies|
    dependencies.flatten.collect do |dep|
      {:inputs => options.merge(:network => dep), :jobname => dep.clean_name }
    end
  end
  task :database_sweep_benchmark => :tsv do
    regulome_tasks = dependencies.select{|d| d.task_name.to_s == 'regulome' }
    bench_tasks = dependencies.select{|d| d.task_name.to_s == 'decouple_R_benchmark' }

    databases = regulome_tasks.collect do |dep|
      database = dep.recursive_inputs[:databases].first
      if database.nil? || database.empty?
        database = "CollecTRI" 
      end
      database
    end

    databases << 'DoRothEA'

    bench_tasks.zip(databases).inject(nil) do |acc,p|
      dep, database = p
      next acc if dep.error?
      tsv = dep.load
      tsv.add_field "Database" do 
        database
      end

      if acc.nil?
        acc = tsv
      else
        offset = acc.length + 1
        tsv.monitor = true
        tsv.through do |id,values|
          acc[(id.to_i + offset).to_s] = values
        end
      end
      acc
    end
  end

  
  dep :regulome, :databases => :regulome, :compute => :produce do |jobname,options|
    databases = SaezLab.ExTRI_databases
    databases.collect do |db|
      {:inputs => options.merge(:databases => databases - [db]), :jobname => "Remove #{db}" }
    end
  end
  dep :decouple_R_benchmark, :network => :sweep, :canfail => true do |jobname,options,dependencies|
    dependencies.flatten.collect do |dep|
      {:inputs => options.merge(:network => dep), :jobname => dep.clean_name }
    end
  end
  task :database_remove_sweep_benchmark => :tsv do
    regulome_tasks = dependencies.select{|d| d.task_name.to_s == 'regulome' }
    bench_tasks = dependencies.select{|d| d.task_name.to_s == 'decouple_R_benchmark' }

    all_databases = SaezLab.ExTRI_databases
    databases = regulome_tasks.collect do |dep|
      database = (all_databases - dep.recursive_inputs[:databases]).first
      if database.nil? || database.empty?
        database = "CollecTRI" 
      end
      database
    end

    offset = 1
    bench_tasks.zip(databases).inject(nil) do |acc,p|
      dep, database = p
      next acc if dep.error?
      tsv = dep.load
      tsv.add_field "Database" do 
        database
      end

      if acc.nil?
        acc = tsv
        offset = acc.size
      else
        tsv.through do |id,values|
          acc[offset.to_s] = values
          offset += 1
        end
      end

      acc
    end
  end

  #{{{ ++++ R_BENCHMARK BOXPLOT +++

  dep :database_sweep_benchmark
  input :core_dbs, :boolean, "Use only core DBs", false
  input :method, :select, "Method result to plot", :consensus, :select_options => %w(mlm ulm wsum consensus)
  input :metric, :select, "Metric to plot", :mcauroc, :select_options => %w(mcauroc mcauprc)
  extension :png
  task :database_sweep_benchmark_boxplot => :binary do |core_dbs,method,metric|
    require 'rbbt/util/R'
    require 'rbbt/util/R/plot'

    tsv = step(:database_sweep_benchmark).load

    boxplot_file = file('boxplot.png')

    tsv = tsv.select("metric" => metric).select("method" => "#{method}_estimate")
    tsv = tsv.select("Database" => %w(CollecTRI ExTRI DoRothEA)) if core_dbs

    R::PNG.ggplot(self.tmp_path, tsv, <<-EOF, 4, 4) 
plot = ggplot(data) + geom_boxplot(aes(x=reorder(Database,score), y=score)) + theme_classic() + rbbt.ggplot2.rotate_x_labels() + xlab('Database') + ylab('#{method} #{metric}')

    EOF
    
    nil
  end

  dep :database_remove_sweep_benchmark
  input :core_dbs, :boolean, "Use only core DBs", false
  input :method, :select, "Method result to plot", :consensus, :select_options => %w(mlm ulm wsum consensus)
  input :metric, :select, "Metric to plot", :mcauroc, :select_options => %w(mcauroc mcauprc)
  extension :png
  task :database_remove_sweep_benchmark_boxplot => :binary do |core_dbs,method,metric|
    require 'rbbt/util/R'
    require 'rbbt/util/R/plot'

    tsv = step(:database_remove_sweep_benchmark).load

    boxplot_file = file('boxplot.png')

    tsv = tsv.select("metric" => metric).select("method" => "#{method}_estimate")
    tsv = tsv.select("Database" => %w(CollecTRI ExTRI DoRothEA)) if core_dbs

    R::PNG.ggplot(self.tmp_path, tsv, <<-EOF, 4, 4) 
plot = ggplot(data) + geom_boxplot(aes(x=reorder(Database,score), y=score)) + theme_classic() + rbbt.ggplot2.rotate_x_labels() + xlab('Removed Database') + ylab('#{method} #{metric}')

    EOF
    
    nil
  end
end
