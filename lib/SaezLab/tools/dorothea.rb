require 'rbbt-util'
module DoRothEA

  Rbbt.claim Rbbt.share.databases.DoRothEA.regulome, :proc do |file|
    require 'rbbt/util/R'
    TmpFile.with_file do |tmp|
      R.run <<-EOF, :monitor => true
rbbt.require('dorothea')
net = dorothea::dorothea_hs
rbbt.tsv.write('#{tmp}', net)
      EOF
      TSV.open(tmp, :type => :list)
    end
  end
end

if __FILE__ == $0
  Log.with_severity 0 do
    Log.tsv Rbbt.share.databases.DoRothEA.regulome.tsv
  end

end
