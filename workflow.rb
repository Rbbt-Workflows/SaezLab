require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/SaezLab'

module SaezLab
  extend Workflow

end

require 'SaezLab/tasks/knocktf.rb'
require 'SaezLab/tasks/decoupler.rb'
require 'SaezLab/tasks/sophia_mueller.rb'
require 'SaezLab/tasks/sweep.rb'

#require 'rbbt/knowledge_base/SaezLab'
#require 'rbbt/entity/SaezLab'

