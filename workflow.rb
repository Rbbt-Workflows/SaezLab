require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/SaezLab'

Workflow.require_workflow "ExTRI"
module SaezLab
  extend Workflow

  input :csv, :file, "CSV file"
  task :csv_to_regulome => :tsv do |csv|
    tsv = TSV.setup({}, "ID~source,target,weight#:type=:list")
    id = 0
    TSV.traverse csv, :type => :array do |line|
      next if line =~ /source/
      id += 1
      source, target, weight, pmids = line.split(",")
      tsv[id] = [source, target, weight]
    end
    tsv
  end

end

require 'SaezLab/tasks/decoupler.rb'
require 'SaezLab/tasks/carnival.rb'
require 'SaezLab/tasks/progeny.rb'
#require 'rbbt/knowledge_base/SaezLab'
#require 'rbbt/entity/SaezLab'

