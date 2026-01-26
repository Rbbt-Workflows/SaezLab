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

  KnockTF.claim KnockTF.expression_new, :proc do
    url = "http://www.licpathway.net/KnockTFv2/public/download_anno/knocktf_v2_main_human.txt"
    url = Rbbt.data['knocktf_v2_main_human.txt']
    data = {}
    datasets = []
    dataset_tf = {}
    fields = nil
    p_dataset = 0
    p_tf = 1
    p_gene = 2
    p_fc = 5
    TSV.traverse url, bar: true, type: :array do |line|
      parts = line.split("\t")
      if fields.nil?
        fields = parts
        next
      end
      dataset = parts[p_dataset]
      tf = parts[p_tf]
      gene = parts[p_gene]
      fc = parts[p_fc]
      dataset = [dataset, tf] * ":"
      datasets << dataset unless datasets.include? dataset
      data[gene] ||= []
      data[gene] << fc
    end

    tsv = TSV.setup({}, key_field: "Associated Gene Name", fields: datasets, type: :list, cast: :to_f )

    data.each do |gene,values|
      tsv[gene] = values
    end
    tsv
  end
end

if __FILE__ == $0
  Log.with_severity 0 do
    Log.tsv KnockTF.expression_new.produce(true).tsv
  end
end
