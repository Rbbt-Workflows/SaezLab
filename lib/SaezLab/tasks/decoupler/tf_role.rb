module SaezLab

  task :tf_role_annotations => :tsv do
    tsv = Rbbt.data.annotations["tf_annotation_activators_repressors.csv"].tsv :sep => ',', :header_hash => '', :type => :list
    tsv.key_field = "Associated Gene Name"

    astrid_activators = Rbbt.data.annotations["Astrid_go_activator_proteins.txt"].list
    astrid_repressors = Rbbt.data.annotations["Astrid_go_repressor_proteins.txt"].list

    tsv.add_field "GO_Activator_0045944_with_children" do |k|
      astrid_activators.include?(k) ? "1" : "0"
    end

    tsv.add_field "GO_Repressor_0000122_with_children" do |k|
      astrid_repressors.include?(k) ? "1" : "0"
    end

    soto = Rbbt.data.annotations["soto2021.csv"].tsv :header_hash => '', :type => :list, :fields => ["Domain type","Activity (H, M or L)", "Confidence (H, M or L)"], :key_field => "TF name"
    soto.key_field = "Associated Gene Name"

    soto.fields = ["Soto Role", "Soto Activity", "Soto Confidence"]

    soto.process "Soto Role" do |v|
      case v
      when "AD"
        "+"
      when "RD"
        "-"
      else
        nil
      end
    end

    taipale = Rbbt.data.annotations["taipale2022.csv"].tsv :header_hash => '', :type => :list
    taipale.key_field = "Associated Gene Name"

    taipale.add_field "Taipale Active Fraction Role" do |k,v|
      s = v.first.to_f
      if s > 1
        "+"
      elsif s < 1
        "-"
      else
        nil
      end
    end

    taipale.add_field "Taipale Intensity Role" do |k,v|
      s = v[1].to_f
      if s > 1
        "+"
      elsif s < 1
        "-"
      else
        nil
      end
    end


    tsv = tsv.attach taipale, :complete => true, :fields => ["Taipale Active Fraction Role", "Taipale Intensity Role"]
    tsv = tsv.attach soto, :complete => true

    krab_super = Rbbt.data.annotations["KRAB_superfam_HS_IPRO3651.txt"].tsv(:header_hash => "").keys
    krab_ancient = Rbbt.data.annotations["KRAB_ancient_HS_IPRO03655.txt"].tsv(:header_hash => "").keys

    uni2name = Organism.identifiers(Organism.default_code("Hsa")).index :target => "Associated Gene Name", :fields => "UniProt/SwissProt Accession", :persist => true

    krab_super_names = uni2name.chunked_values_at krab_super
    krab_ancient_names = uni2name.chunked_values_at krab_ancient

    krab_genes = (krab_super_names + krab_ancient_names).compact.uniq

    krab = TSV.setup(krab_genes, "Associated Gene Name~#:type=:list")

    krab.add_field "KRAB Superfamily" do |k|
      krab_super_names.include?(k)? "true" : nil 
    end

    krab.add_field "KRAB Ancient" do |k|
      krab_ancient_names.include?(k)? "true" : nil 
    end
    Log.tsv krab


    tsv.attach krab, :complete => true
    tsv
  end

  dep :tf_role_annotations
  task :tf_role_manual => :tsv do 
    new = TSV.setup({}, "Associated Gene Name~Sign#:type=:single")
    step(:tf_role_annotations).load.through do |k,values|
      values.collect!{|v| v.empty? ? nil : v}

      if values["KRAB Superfamily"] == "true" && ! values["KRAB Ancient"] == "true"
        sign = "-"
      elsif values["Soto Role"] && values["Soto Confidence"] == "H"
        sign = values["Soto Role"]
      elsif values["Soto Role"] && values["Soto Confidence"] == "M"
        sign = values["Soto Role"]
      elsif values["Taipale Active Fraction Role"]
        sign = values["Taipale Active Fraction Role"]
      elsif values["Soto Role"]
        sign = values["Soto Role"]
      else
        pos_sum = values.values_at(*%w(KW_Activator GO_Activator_0045944_with_children GO_Activator_0001228)).select{|v| v == "1" }.length
        neg_sum = values.values_at(*%w(KW_Repressor GO_Repressor_0000122_with_children GO_Repressor_0001227)).select{|v| v == "1" }.length
        next if pos_sum == 0 && neg_sum == 0
        sign = neg_sum >= pos_sum ? "-" : "+"
      end

      new[k] = sign
    end
    new
  end

  task :tf_role => :tsv do
    tsv = Rbbt.data.annotations.Astrid_TF_role["CollecTRI_TF-role_new_draft_211122.tsv"].tsv :type => :list
    tsv.key_field = "Associated Gene Name"
    tsv.add_field "Sign" do |k,v|
      case v["STRICT_agreement (GO/UniProt-StructureFunction)"]
      when "Act"
        "+"
      when "Repr"
        "-"
      else
        nil
      end
    end

    tsv = tsv.select("Sign" => %w(+ -))
    tsv.slice(["Sign"]).to_single
  end

end
