module SaezLab

  input :expression, :tsv, "Gene expression levels"
  input :organism_name, :select, "Organism name", "Human", :select_options => %w(Human Mouse)
  task :progeny => :tsv do |expression,org|
    require 'rbbt/util/R'

    expression.R <<-EOF
rbbt.require('progeny')
data = subset(data,rowSums(is.na(data))==0)
data = progeny(as.matrix(data), scale=TRUE, organism=#{R.ruby2R org})
    EOF
  end
end
