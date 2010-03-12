require 'rubygems'
require 'soap/wsdlDriver'
require 'base64'

WSDL_URI  = File.join(File.dirname(File.expand_path(__FILE__)), 'wsdl', 'TaLassoWS.wsdl')
DATA_PATH = File.join(File.dirname(File.expand_path(__FILE__)), '../data/examples/Endotelio')
server =  SOAP::WSDLDriverFactory.new(WSDL_URI).create_rpc_driver 

data_id = server.upload(
  File.open(File.join(DATA_PATH, 'geneExpression.txt')).read, 
  File.open(File.join(DATA_PATH, 'mirnaExpression.txt')).read, 
  File.open(File.join(DATA_PATH, 'gene.txt')).read, 
  File.open(File.join(DATA_PATH, 'mirna.txt')).read, 
  File.open(File.join(DATA_PATH, 'samples.txt')).read
) 


job = server.lasso(data_id, %w(tarbase), '')

while ! server.done job
  puts "."
  sleep 5
end

raise "Error: " + server.messages(job).last if server.error job

results = server.results(job)

File.open('targets.txt', 'w') do |f| f.write Base64.decode64 server.result(results[0]) end
File.open('gene.txt', 'w') do |f| f.write Base64.decode64 server.result(results[1]) end
File.open('mirna.txt', 'w') do |f| f.write Base64.decode64 server.result(results[2]) end
