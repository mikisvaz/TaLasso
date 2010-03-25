$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
require 'TaLasso'
require 'validation'

GOLDSTANDARDS = %w(tarbase)

helper :validate do |job_name, goldstandard|
  info[:hypergeometric] ||= {}
  info[:hits]           ||= {}

  step(:validation, "Loading #{goldstandard} for validation")

  results, positions = TaLasso::Validation.positions(File.join(workdir, job_name), File.join(TaLasso::PUTATIVE_DIR, goldstandard)).values_at(:results, :positions)

  step(:validation, "Performing hypergeometric validation with #{goldstandard}")

  pvalues = TaLasso::Validation.validate(positions, results, File.join(workdir,job_name,"#{goldstandard}_pvalues.png"))
  info[:hypergeometric][goldstandard.to_sym] = pvalues
  
  step(:validation, "Performing rank validation with #{goldstandard}")
  
  score   = TaLasso::Validation.hit_score(positions, results, File.join(workdir,job_name,"#{goldstandard}_hits.png"))
  info[:hits][goldstandard.to_sym] = score
end

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
   '{JOB}/mirna.txt',
   '{JOB}/tarbase_pvalues.png',
   '{JOB}/tarbase_hits.png',
] do |data_id, putatives|
  info[:type] = 'lasso'

  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing Lasso optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.lasso(output_path, input_path, putatives)
  
  GOLDSTANDARDS.each do |goldstandard| validate(job_name, goldstandard) end
end

task :correlation, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt',
   '{JOB}/tarbase_pvalues.png',
   '{JOB}/tarbase_hits.png',
] do |data_id, putatives|
  info[:type] = 'correlation'

  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Calculating correlation")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.correlation(output_path, input_path, putatives)
  
  GOLDSTANDARDS.each do |goldstandard| validate(job_name, goldstandard) end
end

task :gen_mir, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt',
   '{JOB}/tarbase_pvalues.png',
   '{JOB}/tarbase_hits.png',
] do |data_id, putatives|
  info[:type] = 'gen_mir'


  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing GenMir optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.gen_mir(output_path, input_path, putatives)
  
  GOLDSTANDARDS.each do |goldstandard| validate(job_name, goldstandard) end
end

task :simple_gen_mir, %w(data_id putatives), {:putatives => :array}, 
  ['{JOB}/targets.txt',
   '{JOB}/gene.txt', 
   '{JOB}/mirna.txt',
   '{JOB}/tarbase_pvalues.png',
   '{JOB}/tarbase_hits.png',
] do |data_id, putatives|
  info[:type] = 'simple_gen_mir'


  output_path = File.join(workdir, job_name)
  input_path  = File.join(workdir, data_id)

  step(:process, "Performing simple GenMir optimization")
  FileUtils.mkdir_p output_path unless File.exists? output_path
  TaLasso.simple_gen_mir(output_path, input_path, putatives)
  
  GOLDSTANDARDS.each do |goldstandard| validate(job_name, goldstandard) end
end
