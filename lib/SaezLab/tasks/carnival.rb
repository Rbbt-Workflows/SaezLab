module SaezLab

  input :knockouts, :array, "Knockout genes"
  input :knockins, :array, "Knockins genes"
  input :signaling_network, :tsv, "Signaling network"
  input :tf_activities, :tsv, "Transcription factor activity values"
  task :carnival => :tsv do |knockouts,knockins,network,tf_activities|
    require 'rbbt/util/R'

    work = file('work')
    net = work.net
    tf = work.tf
    prog = work.progeny

    network_sif = file('network.sif')
    activity = file('activity.tsv')


    Open.write(net, network.to_list.to_s)
    Open.write(tf, tf_activities.to_s)

    Misc.in_dir file('work') do
      R.run <<-EOF
rbbt.require('CARNIVAL')

knockouts = #{R.ruby2R(knockouts || [])}
knockins = #{R.ruby2R(knockins || [])}
perturbations = c(rep(-1, length(knockouts)), rep(1, length(knockins)))
names(perturbations) = c(knockouts, knockins)

net_tmp = rbbt.tsv('#{net}', row.names=NULL)
net = data.frame(source=net_tmp$source, interaction=as.numeric(net_tmp$interaction), target=net_tmp$target)

tf = rbbt.tsv('#{tf}')
meas = as.matrix(tf)[,1]

res = runCARNIVAL(inputObj = perturbations, measObj=meas, netObj=net)

rbbt.tsv.write('#{network_sif}', res$weightedSIF)
rbbt.tsv.write('#{activity}', res$attributesAll[[1]])
      EOF
    end

    activity.tsv :key_field => "Nodes", :fields => ["Activity"], :type => :single, :cast => :to_f
  end
end
