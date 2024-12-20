require 'SaezLab/sources/knocktf'
module SaezLab
  task :knocktf_old => :tsv do

    tsv = nil
    genes = nil
    gene_values = nil
    dataset = nil
    TSV.traverse Rbbt.data.knockTF["knockTF_raw_table.txt"], :type => :array, :bar => self.progress_bar do |line|
      next if line =~ /^Sample_ID/
      sample, tf, gene, treat, control, fc, log2fc, rank, pvalue, up_down = line.split(/\t/)

      if dataset != sample
        if gene_values
          if tsv.nil?
            genes = gene_values.keys.sort
            tsv = TSV.setup({}, :key_field => "Dataset", :fields => genes, :type => :list, :cast => :to_f)
          end
          tsv[dataset] = gene_values.chunked_values_at(genes).collect{|v| v.nil? ? 0 : v }
        end
        dataset = sample
        gene_values = {}
      end

      log2fc = 0 if log2fc.empty? || log2fc.nil? ||  log2fc == ""
      
      gene_values[gene] = log2fc 
    end

    tsv
  end

  task :knocktf_ground_truth_old => :tsv do
    Rbbt.data.knockTF["knockTF_raw_table.txt"].tsv :fields => ["TF"], :header_hash => "", :type => :single
  end

  #input :knocktf_logfc_threshold, :float, "Threshold to select valid experiments", -1
  task :knocktf => :tsv do |threshold|
    #metadata = KnockTF.metadata.tsv 
    #good = metadata.select("logFC"){|lfc| lfc.to_f <= threshold }.keys

    #log :tfs, metadata.column("TF").values_at(*good).uniq.length

    expression = KnockTF.expression.tsv
    #expression.slice(good).transpose("Sample_ID")
    expression.transpose("Sample_ID")
  end

  input :knocktf_logfc_threshold, :float, "Threshold to select valid experiments", -1
  task :knocktf_ground_truth => :tsv do |threshold|
    metadata = KnockTF.metadata.tsv 
    good = metadata.select("logFC"){|lfc| lfc.to_f <= threshold }.keys
    metadata.column("TF").subset(good).to_single
  end

  input :dorothea_confidence, :select, "Confidence level for DoRothEA", "C", :select_options => %w(A B C D E)
  task :dorothea => :tsv do |confidence|
    tsv = Rbbt.share.databases.DoRothEA.regulome.tsv :type => :list
    tsv = tsv.select "confidence" do |c| c <= confidence end
    tsv.fields = %w(source confidence target weight)
    tsv
  end

  dep :knocktf, :jobname => "Default"
  input :dorothea, :boolean, "Use DoRothEA instead of CombTRI literature counts regulome", false
  input :dorothea_confidence, :select, "Confidence level for DoRothEA", "C", :select_options => %w(A B C D E)
  dep :manual_regulome do |jobname,options|
    if options[:manual_regulome]
      {:task => :manual_regulome, :inputs => options}
    else
      nil
    end
  end
  dep :dorothea do |jobname,options|
    if options[:dorothea]
      {:task => :dorothea, :inputs => options}
    else
      nil
    end
  end
  dep :regulome do |jobname,options|
    if options[:dorothea]
      dependencies.select{|d| d.task_name == :dorothea }.first
    elsif options[:manual_regulome]
      dependencies.select{|d| d.task_name == :manual_regulome }.first
    else
      {:task => :regulome, :inputs => options}
    end
  end
  dep_task :decoupler_knocktf, SaezLab, :decoupler, :matrix => :knocktf, :network => :regulome do |jobname,options,dependencies|
    regulome = dependencies.flatten.last
    {:inputs => options.merge(:network => regulome) }
  end

  dep :decoupler_knocktf
  dep_task :decoupler_knocktf_status, SaezLab, :decoupler_status, "SaezLab#decoupler" => :decoupler_knocktf

  dep :decoupler_knocktf
  dep_task :decoupler_knocktf_congruent, SaezLab, :decoupler_congruent_tagets, "SaezLab#decoupler" => :decoupler_knocktf

  dep :decoupler_knocktf
  dep_task :decoupler_knocktf_incongruent, SaezLab, :decoupler_incongruent_tagets, "SaezLab#decoupler" => :decoupler_knocktf

  dep :decoupler_knocktf
  dep :regulome
  dep_task :decoupler_knocktf_congruency, SaezLab, :decoupler_congruency, "SaezLab#decoupler" => :decoupler_knocktf


  dep :decoupler_knocktf
  dep :knocktf_ground_truth, :jobname => "Default"
  task :decoupler_predictions_knocktf => :tsv do
    ground_truth = step(:knocktf_ground_truth).load
    result = step(:decoupler_knocktf).load

    fields = result.fields

    res = TSV.setup({}, "Dataset~TF,Value,Rank,pvalue#:type=:list")
    ground_truth.each do |dataset,tf|
      dataset_values = result[dataset]
      begin
        if fields.include? "#{tf} (pvalue)"
          pvalue = dataset_values["#{tf} (pvalue)"]
        else
          pvalue = nil
        end

        if fields.include? "#{tf} (consensus_estimate)"
          value = dataset_values["#{tf} (consensus_estimate)"]

          all_values = dataset_values.fields.zip(dataset_values).
            select{|f,v| f.include? "(consensus_estimate)" }.
            collect{|f,v| v ? v.to_f : nil }.compact.sort
        else
          value = dataset_values[tf]
          all_values = dataset_values.fields.zip(dataset_values).
            select{|f,v| ! f.include? "(" }.
            collect{|f,v| v ? v.to_f : nil }.compact.sort
        end
      rescue
        Log.warn $!.message
        value = nil
      end
      rank = all_values.index(value).to_f / all_values.length if value
      res[dataset] = [tf, value, rank,pvalue]
    end

    res
  end

  dep :decoupler_predictions_knocktf
  input :activation_threshold, :float, "Threshold for activation", 0
  input :pvalue_threshold, :float, "Threshold for activation", 0.1
  task :evaluate_knocktf => :tsv do |activation_threshold,pvalue_threshold|
    tsv = TSV.setup({}, "Statistic~Value#:type=:single")
    tsv["Missing"]      ||= 0
    tsv["Miss-matches"] ||= 0
    tsv["Matches"]      ||= 0

    ranks = []
    TSV.traverse step(:decoupler_predictions_knocktf) do |dataset, values|
      tf, value, rank, pvalue = values
      ranks << rank.to_f if rank || ! rank.empty?

      value = nil if value.to_f.abs < activation_threshold
      value = nil if pvalue.to_f > pvalue_threshold

      tsv["Missing"]      += 1 if value.nil? || value.empty?
      tsv["Miss-matches"] += 1 if value && value.to_f > 0
      tsv["Matches"]      += 1 if value && value.to_f < 0
    end

    tsv["Accuracy %"] = (100.0 * tsv["Matches"] / (tsv["Miss-matches"] + tsv["Matches"])).round(2)
    tsv["Average rank"] = Misc.mean(ranks).round(4)

    tsv
  end

  dep :regulome
  dep_task :evaluate_knocktf_auroc, SaezLab, :decouple_R_benchmark_auroc, :network => :regulome

  dep :decoupler_knocktf
  task :activity_expression => :tsv do
    decoupler = step(:decoupler_knocktf).load
    log2 = step(:knocktf).load
    genes = decoupler.fields.collect{|f| f.sub(/ \(p.*\)/,'') }.uniq
    genes &= log2.fields

    data = TSV::Dumper.new(:key_field => "ID", :fields => "Associated Gene Name,Activity,Expression,Dataset,pvalue".split(","), :type => :list)
    data.init

    id = 0
    TSV.traverse decoupler, :into => data, :bar => self.progress_bar do |dataset,decoupler_values|
      log2_values = log2[dataset]
      next if log2_values.nil?

      res = []
      res.extend MultipleResult
      genes.each do |gene|
        activity = decoupler_values[gene]
        activity_pvalue = decoupler_values[gene + " (pvalue)"]
        expression = log2_values[gene]
        next if activity_pvalue > 0.1
        next if expression == 0
        id += 1
        res << [id, [gene, activity, expression, dataset,activity_pvalue]]
      end
      res
    end
  end

  dep :activity_expression
  task :activity_expression_consistency => :tsv do
    genes = TSV.setup({}, "Associated Gene Name~Matches,Miss-matches#:type=:list#:cast=:to_i")

    TSV.traverse step(:activity_expression) do |id,values|
      gene, activity, expression, dataset, pvalue = values

      genes[gene] ||= [0, 0]
      if activity.to_f > 0 == expression.to_f > 0
        genes[gene][0] += 1 
      else
        genes[gene][1] += 1
      end
    end
    genes
  end

  dep :activity_expression_consistency
  task :activity_expression_consistency_ratio => :float do
    data = step(:activity_expression_consistency).load 
    matches = Misc.sum(data.column("Matches").values)
    miss_matches = Misc.sum(data.column("Miss-matches").values)
    matches / miss_matches
  end


end
