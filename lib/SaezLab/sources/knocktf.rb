require 'rbbt-util'
require 'rbbt/resource'
require 'rbbt/tsv/csv'
require 'rbbt/sources/organism'

module KnockTF
  extend Resource

  self.subdir = 'share/databases/KnockTF'

  def self.organism
    Organism.default_code("Hsa")
  end

  KnockTF.claim KnockTF.expression, :proc do
    tsv = TSV.csv("https://zenodo.org/record/7035528/files/knockTF_expr.csv?download=1", cast: :to_f)
    tsv.key_field = "Dataset"
    tsv.namespace = KnockTF.organism
    tsv.transpose("Associated Gene Name")
  end

  KnockTF.claim KnockTF.metadata, :proc do
    TSV.csv("https://zenodo.org/record/7035528/files/knockTF_meta.csv?download=1")
  end
end

if __FILE__ == $0
  Log.with_severity 0 do
    Log.tsv KnockTF.expression.tsv
    Log.tsv KnockTF.metadata.tsv
  end
end
