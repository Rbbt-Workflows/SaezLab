module SaezLab

  dep :evaluate_knocktf, :negative_evidence_proportion => :sweep,
    :sign_support_balance => :sweep, :strict_negative => :sweep, :support_evidence_max => :sweep do |jobname,options|

    sweep = [0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1]
    support_evidence_max_sweep = [3,5,10]
    negative_evidence_proportion = 1

    jobs = []
    sweep.each do |sign_support_balance|
      [0.1,0.2,0.3,0.5].each do |strict_negative|
        support_evidence_max_sweep.each do |support_evidence_max|
          job_options = options.merge(
            :negative_evidence_proportion => negative_evidence_proportion,
            :sign_support_balance => sign_support_balance,
            :strict_negative => strict_negative,
            :support_evidence_max => support_evidence_max)

          jobs << {:inputs => job_options}
          #jobs << SaezLab.job(:evaluate_knocktf, nil, job_options)
        end
      end
    end

    jobs
  end
  task :sweep => :tsv do
    fields = [:negative_evidence_proportion, :sign_support_balance, :strict_negative, :support_evidence_max]
    dependencies.inject(nil) do |acc,dep|
      expression = fields.zip(dep.recursive_inputs.values_at *fields).collect{|f,v| [f,v] * "="} * ","
      tsv = dep.load.transpose("Experiment")
      tsv[expression] = tsv.delete "Value"
      
      acc = acc.nil? ? tsv : acc.merge!(tsv)
    end
  end

  dep :evaluate_knocktf, :databases => :placeholder do |jobname,options|
    Workflow.require_workflow "ExTRI"
    extri = ExTRI.job('pairs_final').load
    databases = extri.fields.collect do |f,v|
      m = f.match(/\[(.*)\] (.*)/)
      next unless m
      db, field = m.values_at 1, 2
      db
    end.flatten.uniq

    databases.collect do |db|
      {:inputs => options.merge(:databases => [db]) }
    end
  end
  dep :evaluate_knocktf, :databases => nil, :dorothea => true
  task :database_sweep => :tsv do
    dependencies.inject(nil) do |acc,dep|
      if dep.recursive_inputs[:dorothea]
        database = "Dorothea"
      else
        database = dep.recursive_inputs[:databases].first
        if database.nil? || database.empty?
          database = "All" 
        end
      end
      tsv = dep.load.transpose("Database")
      tsv[database] = tsv.delete "Value"
      
      acc = acc.nil? ? tsv : acc.merge!(tsv)
    end
  end
end
