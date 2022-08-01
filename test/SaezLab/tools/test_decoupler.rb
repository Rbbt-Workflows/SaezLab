require File.join(File.expand_path(File.dirname(__FILE__)), '../..', 'test_helper.rb')
require 'SaezLab/tools/decoupler'

class TestDecoupler < Test::Unit::TestCase
  def _test_decoupler
    mat = net = nil
    SaezLab._pycall_decoupler do
      mat, net = decoupler.get_toy_data()
    end

    acts, acts_normal, pvals = SaezLab.decoupler(:gsea, mat, net, min_n: 0, times: 100)
    assert acts.keys.include?("S01")

    acts, acts_normal, pvals = SaezLab.decoupler(:viper, mat, net, min_n: 0)
    assert acts.keys.include?("S01")
    assert_nil pvals
  end

  def test_tsv
    mat = net = nil
    SaezLab._pycall_decoupler do
      mat, net = decoupler.get_toy_data()
    end

    mat_tsv = TSV.open(StringIO.new(RbbtPython.df2tsv(mat).to_s))
    Log.tsv mat_tsv
    mat = RbbtPython.tsv2df(mat_tsv)

    acts, acts_normal, pvals = SaezLab.decoupler(:gsea, mat, net, min_n: 0, times: 100)
    assert acts.keys.include?("S01")
  end
end

