require 'rbbt/util/R'

module SaezLab

  task :resource_aggregation => :array do
    workdir = file('work')
    Misc.in_dir workdir do 
      Open.mkdir workdir.data
      Open.mkdir workdir.output
      R.run Rbbt.modules.CollecTRI.scripts.CollecTRI["01_resource_aggregation.R"].read
    end
    workdir.glob("*")
  end

  dep :resource_aggregation
  task :mode_of_regulation => :tsv do
  end

end
