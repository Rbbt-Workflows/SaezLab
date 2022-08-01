require 'SaezLab/tools/decoupler'

module SaezLab

  input :method, :select, "Method to run", :ulm, :select_options => %w(ulm mln wsum combined)
  input :matrix, :tsv, "Matrix of readouts"
  input :network, :tsv, "Network with weights"
  input :times, :integer, "Not always used",nil
  input :min_n, :integer, "Minimum of target features per biological entity ", 5
  task :decoupler => :tsv do |method,matrix,network,times,min_n|

    options = {times: times, min_n: min_n}
    options.delete_if{|k,v| v.nil? }
    acts, norm_acts, pvals = SaezLab.decoupler(method, matrix, network, options)

    Open.write(file('acts'), acts.to_s) if acts
    Open.write(file('norm_acts'), norm_acts.to_s) if norm_acts
    Open.write(file('pvals'), pvals.to_s) if pvals

    norm_acts.fields = norm_acts.fields.collect{|f| "#{f} (norm)" } if norm_acts
    pvals.fields = pvals.fields.collect{|f| "#{f} (pvalue)" } if pvals

    acts = acts.attach norm_acts if norm_acts
    acts = acts.attach pvals if pvals

    acts
  end

  dep :decoupler
  task :decoupler_status => :array do 

    matrix = self.recursive_inputs[:matrix]
    matrix = TSV.open matrix unless TSV === matrix
    matrix.unnamed = true

    activations = step(:decoupler).load
    good_fields = activations.fields.select{|f| ! f.include?("(") || f.include?("(combined)") }
    activations = activations.slice good_fields
    activations.fields = good_fields.collect{|f| f.split(" ").first}
    activations.unnamed = true

    network = self.recursive_inputs[:network]
    network = TSV.open network unless TSV === network
    network = network.to_double

    network = Association.open network, :source => "source", :target => "target", :merge => true, :fields => %w(weight), :persist => true
    network.unnamed = true

    associations = Association.index network
    associations.unnamed = true

    all_sources = network.keys

    output = file('output')

    cpus = config :cpus, :decoupler, :default => 3

    sources = activations.fields
    TSV.traverse activations, :bar => self.progress_bar("processing samples"), :cpus => cpus do |sample,activations|
      target_source_status = TSV.setup({}, :key_field => "target", :fields => sources, :type => :list)

      source_activations = Misc.zip2hash(sources, activations)
      i = 0
      TSV.traverse sources do |source|
        activation = source_activations[source]
        targets = network[source].first

        targets.each do |target|
          pair = [source, target] * "~"

          next unless associations.include?(pair)

          weight = associations[pair].first.to_f

          score = weight.to_f * activation.to_f

          if score.abs < 0.1
            status = "neutral"
          else
            status = score > 0 ? "match" : "missmatch"
          end

          target_source_status[target] ||= [nil] * sources.length
          target_source_status[target][i] = status
        end
        i += 1
      end

      Open.write output[sample], target_source_status.to_s
    end

    output.glob("*")
  end

  dep :decoupler_status
  task :decoupler_congruent_tagets => :tsv do
    files = step(:decoupler_status).file('output').glob("*")
    datasets = files.collect{|f| File.basename(f) }.sort

    congruent = TSV.setup({}, :key_field => "target", :fields => datasets, :type => :double)

    congruent.unnamed = true

    i = 0
    TSV.traverse files, :bar => self.progress_bar("Processing files") do |file|
      TSV.traverse file do |target,status,sources|
        congruent[target] ||= [[]] * datasets.length

        target_dataset_matches = congruent[target]
        current = target_dataset_matches[i].dup
        Misc.zip_fields([sources, status]).each do |source, stat|
          if stat == "match"
            current << source
          end
        end
        target_dataset_matches[i] = current
        congruent[target] = target_dataset_matches
      end
      i += 1
    end

    congruent
  end

  dep :decoupler_status
  task :decoupler_incongruent_tagets => :tsv do
    files = step(:decoupler_status).file('output').glob("*")
    datasets = files.collect{|f| File.basename(f) }.sort

    congruent = TSV.setup({}, :key_field => "target", :fields => datasets, :type => :double)

    congruent.unnamed = true

    i = 0
    TSV.traverse files, :bar => self.progress_bar("Processing files") do |file|
      TSV.traverse file do |target,status,sources|
        congruent[target] ||= [[]] * datasets.length

        target_dataset_matches = congruent[target]
        current = target_dataset_matches[i].dup
        Misc.zip_fields([sources, status]).each do |source, stat|
          if stat == "missmatch"
            current << source
          end
        end
        target_dataset_matches[i] = current
        congruent[target] = target_dataset_matches
      end
      i += 1
    end

    congruent
  end

  dep :decoupler_status
  dep :decoupler_congruent_tagets
  dep :decoupler_incongruent_tagets
  task :decoupler_congruency => :tsv do

    congruent   = step(:decoupler_congruent_tagets).load
    incongruent = step(:decoupler_incongruent_tagets).load
    datasets    = congruent.fields

    Workflow.require_workflow "ExTRI"
    extri_pre = ExTRI.job('pairs_final').load
    extri = extri_pre.annotate({})
    extri_pre.each do |k,v| 
      extri[k.sub(':', '~') ] = v
    end

    log :congruent, "Loading congruent pairs"
    congruent_pairs = {}

    congruent.each do |target,dataset_info|
      dataset_info.each_with_index do |list,i|
        dataset = datasets[i]
        list.each do |source|
          pair = [source, target] * "~"
          congruent_pairs[pair] ||= []
          congruent_pairs[pair] << dataset
        end
      end
    end

    log :incongruent, "Processing incongruent into explained and unexplained"
    explained_pairs = {}
    unexplained_pairs = {}

    incongruent.each do |target,dataset_info|
      dataset_info.each_with_index do |list,i|
        dataset = datasets[i]
        list.each do |source|
          pair = [source, target] * "~"
          explanations = congruent[target][i]
          if explanations.any?
            explained_pairs[pair] ||= []
            explained_pairs[pair] << explanations*";"
          else
            unexplained_pairs[pair] ||= []
            unexplained_pairs[pair] << dataset
          end
        end
      end
    end

    log :merging_info, "Attaching a single TSV file"
    TSV.setup(congruent_pairs, :type => :flat, :key_field => "Pair",
              :fields => %w(congruent))
    TSV.setup(unexplained_pairs, :type => :flat, :key_field => "Pair",
              :fields => %w(unexplained))
    TSV.setup(explained_pairs, :type => :flat, :key_field => "Pair",
              :fields => %w(explained))

    tsv =    congruent_pairs.to_double.
      attach(unexplained_pairs.to_double, :complete => true).
      attach(explained_pairs.to_double, :complete => true)

    log :counts, "Turning into counts"

    tsv = tsv.process "congruent" do |list|
      [list.length]
    end

    tsv = tsv.process "unexplained" do |list|
      [list.length]
    end

    tsv = tsv.add_field "explanations" do |v,values|
      list = values["explained"]
      Misc.counts(list.collect{|e| e.split(";") }.flatten).
        sort_by{|tf,c| -c }.
        collect{|tf,c| [tf,c] * "="}
    end

    tsv = tsv.process "explained" do |list|
      [list.length]
    end

    log :regulome, "Attaching regulome"
    network = self.recursive_inputs[:network]
    network = TSV.open network unless TSV === network
    network.unnamed = true

    network = Association.open network.to_double, :merge => true, :persist => true,
      :source => "source", :target => "target",
      :fields => %w(weight)

    network.unnamed = true

    associations = Association.index network
    associations.unnamed = true

    associations.key_field = tsv.key_field
    tsv = tsv.attach associations

    log :ExTRI, "Attaching ExTRI information"
    extri.key_field = tsv.key_field
    extri.unnamed = true

    tsv.attach extri
  end

end
