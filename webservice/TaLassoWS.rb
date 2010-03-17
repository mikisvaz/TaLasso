$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
require 'TaLasso'

serve :upload, %w(geneExpression mirnaExpression genes mirna samples) do |ge, me, gene, mirna, samples|
  name = SimpleWS::Jobs::Scheduler.random_name('data-')
  path = File.join(workdir, name)

  FileUtils.mkdir_p path
  File.open(File.join(path, 'geneExpression.txt'), 'w')  do |f| f.write ge end
  File.open(File.join(path, 'mirnaExpression.txt'), 'w') do |f| f.write me  end
  File.open(File.join(path, 'gene.txt'), 'w')            do |f| f.write gene end
  File.open(File.join(path, 'mirna.txt'), 'w')           do |f| f.write mirna end
  File.open(File.join(path, 'samples.txt'), 'w')         do |f| f.write samples end

  name
end

task :lasso, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt'] do |data_id, putatives|

  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing Lasso optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.lasso(output_path, input_path, putatives)
end

task :correlation, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt'] do |data_id, putatives|

  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing Lasso optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.correlation(output_path, input_path, putatives)
end

task :gen_mir, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt'] do |data_id, putatives|

  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing Lasso optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.gen_mir(output_path, input_path, putatives)
end

task :simple_gen_mir, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt'] do |data_id, putatives|

  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing Lasso optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.simple_gen_mir(output_path, input_path, putatives)
end

