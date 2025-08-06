module SaezLab

  dep ExTRI, :CollecTRI, :jobname => "Default"
  input :high_confidence, :boolean, "Use only high confidence entries", false
  input :databases, :array, "Limit regulome to certain databases", []
  task :count_signs => :tsv do |hc,databases|

    Workflow.require_workflow "ExTRI"

    collecTRI = step(:CollecTRI)

    pair_sign_info = {}
    TSV.traverse collecTRI, :bar => self.progress_bar("Processing pairs file") do |pair,values,fields|
      pair = pair.first if Array === pair
      db_fields = {}

      fields.zip(values).each do |f,v|
        m = f.match(/\[(.*)\] (.*)/)
        next unless m
        db, field = m.values_at 1, 2

        v ||= ""

        field = "Sign" if %w(regulation sign effect activation/repression ).include? field.downcase
        field = "Confidence" if %w(confidence).include? field.downcase
        db_fields[db] ||= {}
        db_fields[db][field] = v
      end

      db_sign_info = {}
      db_fields.each do |db,fields|
        next if databases && databases.any? && ! databases.include?(db)
        present = fields["present"]
        next if present.nil?  || present.empty?

        pmids = fields["PMID"] || []
        signs = fields["Sign"] || []
        confs = fields["Confidence"] || []

        if pmids.empty? && signs.empty?
          pmids = [""]
          signs = ["unknown"]
        end

        pmids = [""] * signs.length if pmids.empty?
        signs = ["unknown"] * pmids.length if signs.empty?
        confs = [""] * pmids.length if confs.empty?

        signs = signs.collect do |s| 
          case s.downcase
          when 'activate', 'up', 'activation', 'positive', 'go:2000144', '+', 'stimulate'
            '+'
          when 'deactivate', 'repress', 'repression', 'down', 'negative', 'go:0033234', '-', 'inhibit'
            '-'
          when '+_-', '-_+'
            '+_-'
          when 'unknown', '', 'NA', 'not_applicable', '?', 'undefined'
            '~'
          else
            raise "What sign: '#{s}'"
          end
        end

        db_sign_info[db] = signs.zip(pmids, confs)
      end

      pair_sign_info[pair] = db_sign_info unless db_sign_info.empty?
    end

    dbs = pair_sign_info.values.first.keys.sort

    tsv = TSV.setup({}, :key_field => "TF-TG", :fields => %w(Sign Evidence), :type => :double)
    pair_sign_info.each do |pair,sign_info|
      sign_evidence = {}
      sign_info.each do |db,info|
        info.each do |s,e,c|
          e = "rand-#{Misc.digest(rand(10000).to_s)[0..4]}" if e.nil? || e.empty?
          next if hc && c && c.downcase == 'low'
          sign_evidence[s] ||= []
          sign_evidence[s] << e.split(";").uniq
        end
      end
      tsv[pair] = Misc.zip_fields(sign_evidence.collect{|s,e| [s, e*";"] })
    end

    tsv
  end

  dep :count_signs
  dep :tf_role, :jobname => "Default"
  task :missing_sign => :array do
    tf_role = step(:tf_role).load

    TSV.traverse step(:count_signs), :into => :stream, :bar => self.progress_bar("Generating regulome") do |pair,values|
      pair = pair.first if Array === pair
      source, target = pair.split(":")

      ### This loop is run for each pair

      # sign_evidence_pmids will hold the list of PMIDS support each sign: +,
      # -, and ~ (unknown or not specified)
      sign_evidence_pmids = {"+" => [], "-" => [], "~" => []}
      Misc.zip_fields(values.reverse).each do |e,s|
        s = "" if s.nil?
        sign_evidence_pmids[s] = e.split(";")
      end

      if sign_evidence_pmids["+"].any? || sign_evidence_pmids["-"].any?
        pos, neg = sign_evidence_pmids.values_at("+", "-").collect{|v| v.length}

        min = [pos, neg].min
        sum = pos + neg
        
        next if min < sum.to_f * 0.3
      end

      next if %w(- +).include? tf_role[source]

      pair
    end
  end

  dep :count_signs
  task mainly_repressor: :array do
    tf_signs = {}
    TSV.traverse step(:count_signs), :bar => self.progress_bar("Generating regulome") do |pair,values|
      source, target = pair.split(":")
      signs, pmids = values
      counts = pmids.collect{|v| v.split(';').uniq.length }
      sign_counts = Misc.zip2hash(signs, counts)
      negative = sign_counts["-"] || 0
      positive = sign_counts["+"] || 0

      tf_signs[source] ||= 0
      if negative >= positive
        tf_signs[source] -= 1
      else 
        tf_signs[source] += 1
      end
    end
    tf_signs.select{|tf,c| c < 0 }.collect{|tf,c| tf }
  end

  dep :count_signs
  dep :mainly_repressor
  dep :tf_role, :jobname => "Default"
  input :negative_evidence_proportion, :float, "(DEACTIVATED) Proportion of negative evidence for a negative sign", 1
  input :support_evidence_max, :integer, "Number of supporting evidence papers maximum to consider", 100
  input :sign_support_evidence_max, :integer, "(DEACTIVATED) Number of supporting evidence papers maximum to consider", 3
  input :sign_support_balance, :float, "Poportion of final weight that comes from sign", 0.1
  input :strict_negative, :float, "Proportion of negative evidence required to make pair negative", 0.5
  input :use_tf_role, :boolean, "Use tf_role to fix negative values", false
  input :auto_regulation_weight, :float, "Set weight for auto-regulation", 0.5
  task :regulome => :tsv do |negative_evidence_proportion,support_evidence_max,sign_support_evidence_max,sign_support_balance,strict_negative,use_tf_role,auto_regulation_weight|


    auto_regulation_weight = nil if auto_regulation_weight.to_f == 0

    tf_role = step(:tf_role).load if use_tf_role
    mainly_repressor = step(:mainly_repressor).load

    dumper = TSV::Dumper.new :key_field => "ID", :fields => %w(source target weight), :type => :list
    dumper.init
    id = 1
    controversies = 0
    begin
      TSV.traverse step(:count_signs), :into => dumper, :bar => self.progress_bar("Generating regulome") do |pair,values|
        pair = pair.first if Array === pair
        source, target = pair.split(":")
        if auto_regulation_weight && source == target
          weight = auto_regulation_weight
        else

          ### This loop is run for each pair

          # sign_evidence_pmids will hold the list of PMIDS support each sign: +,
          # -, and ~ (unknown or not specified)
          sign_evidence_pmids = {"+" => [], "-" => [], "~" => []}
          Misc.zip_fields(values.reverse).each do |e,s|
            s = "" if s.nil?
            sign_evidence_pmids[s] = e.split(";").uniq
          end

          # NOTE: The lists of PMIDs for each sign might include repetitions, so
          # if two databases reviewed the same PMID and reported it, it will be
          # counted twice. Perhaps we can change that

          # Here we just count the pairs in which we have the same PMID used as
          # evidence for contradicting signs (the & sign means intersecting two
          # lists). Not used in scoring.
          if (sign_evidence_pmids["-"] & sign_evidence_pmids["+"]).any?
            controversies += 1
          end

          # Here I take out all PMIDs that indicate some sign from the list of ~
          # because maybe two databases report the same but one didn't care to
          # indicate the sign. In that case I assume that the PMID does indicate
          # the sign and ignore the entry that doesn't reflect that
          sign_evidence_pmids["~"] = sign_evidence_pmids["~"] - (sign_evidence_pmids["+"] + sign_evidence_pmids["-"])

          # Here we turn the lists of PMIDs into counts, which is basically just
          # the length of the list in each sign
          sign_evidence = Hash.new(0)
          sign_evidence_pmids.each do |s,l|
            sign_evidence[s] = l.length
          end

          # Just use some placeholder variables here for clarity. They are turn
          # to continuous variables (float; .to_f) to avoid rounding up when I do operations
          # later with them
          positive_evidence = sign_evidence["+"].to_f
          negative_evidence = sign_evidence["-"].to_f

          signed_evidence = Misc.sum(sign_evidence.values_at("+", "-")).to_f

          non_negative_evidence = Misc.sum(sign_evidence.values_at("+", "~")).to_f
          non_positive_evidence = Misc.sum(sign_evidence.values_at("-", "~")).to_f

          # This is the total number of articles
          total_evidence = Misc.sum(sign_evidence.values).to_f

          # Here we decide if the negative evidence is sufficient to deem the
          # interaction negative: more negative articles than positive, while
          # also the negative evidence divided by positive + unkown must be more than
          # <negative_evidence_proportion>, one of the parameters. Come to think
          # of it the denominator should probably be total_evidence and not
          # non_positive_evidence, otherwise it's not really a proportion. I
          # haven't change this yet
          negative = negative_evidence > positive_evidence && negative_evidence / non_positive_evidence > negative_evidence_proportion

          # As a final rule, the pair is negative if the proportion of negative
          # evidence divided by all the signed_evidence (i.e. ignoring the
          # unknown) is larger than <strict_negative>, another parameter
          negative = true if (negative_evidence / signed_evidence) > strict_negative

          negative = true if mainly_repressor.include?(source) && positive_evidence == 0

          flip = false
          if use_tf_role && tf_role[source] == "-" && positive_evidence == 0
            #iii [pair, sign_evidence_pmids] unless negative
            flip = true unless negative
            negative = true
          end

          # Here we define the sign_weight as the proportion of evidence
          # supporting the chosen sign, which is positive unless it we deemed
          # negative. This means that a totally unkown signed pair will be
          # positive by default, but the sign_weight will be 0
          if negative
            sign = "-"
          else
            sign = "+"
          end

          # Here we compute the support_weight which measure how close we get to
          # having <support_evidence_max> (another parameter) PMIDs not
          # contradicting our sign. So if support_evidence_max is 10 and we have
          # 2 positive, 1 negative and 2 unkown evidences than that would be a
          # positive relationship with 4 non_contradictory_evidences, and the
          # support_weight will be 4 / 10

          non_contradictory_evidence = Misc.sum(sign_evidence.values_at(*[sign, "~"].uniq))
          support_weight = [1, non_contradictory_evidence / support_evidence_max].min

          # ToDo: test if sign_weight now
          correct_sign_evidence = sign_evidence[sign]
          #signed_support_weight = [1, correct_sign_evidence / sign_support_evidence_max].min
          signed_support_weight = [1, correct_sign_evidence / total_evidence].min

          # Here we calculate the final weight by combining sign_weight and
          # support_weight mixed in according to the variable
          # <sign_support_balance>
          weight = (signed_support_weight * sign_support_balance) + (support_weight * (1.0-sign_support_balance))

          # The final weight is made negative when the sign is negative because
          # this is how it is encoded in decopler
          weight = - weight if sign == "-" 

        next if weight.nan?

      end
        id += 1
        [id, [source, target, weight]]
      end
    ensure
      set_info :controversies, controversies
    end
  end

  input :manual_regulome, :file, "Manual regulome" #, Rbbt.data.regulomes["CollecTRI.1-3-23.tsv"].tsv
  task :manual_regulome => :tsv do |manual_regulome|
    TSV === manual_regulome ? TSV.open(manual_regulome) : manual_regulome
  end

end
