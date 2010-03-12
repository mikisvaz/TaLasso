require 'open4'
module TaLasso
  MATLAB='/usr/local/bin/matlab -nojvm -nodesktop -nodisplay'
  MATLAB_LIB= File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'Matlab')
  PUTATIVE_DIR = File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'data', 'targets')

  def self.lasso(output_dir, data_dir, putative, lamb = 1/5)
    paths = putative.collect {|name| File.join(PUTATIVE_DIR, name) }

    FileUtils.mkdir_p output_dir unless File.exists? output_dir
    script =<<-EOF
      addpath('#{MATLAB_LIB}');
      addpath('#{File.join(MATLAB_LIB, 'miLasso')}');
      addpath('#{File.join(MATLAB_LIB, 'miLasso', 'l1_ls_matlab')}');

      [GE, ME, gene, mirna] = loadExpressionData('#{data_dir}');

      paths=[#{paths.collect {|path| "cellstr('" + path + "')" } * ', ' }];

      [GM, geneUnion, mirnaUnion] = loadPutativeData(paths);
      [GM, GE, ME, gene, mirna]   = GMdataIntersect(GM, geneUnion, mirnaUnion, gene, mirna, GE, ME);

      GE = GE - median(GE, 2) * ones(1, size(GE, 2));

      solucion = l1_reg_model(GE, ME, GM, #{ lamb });

      fid = fopen('#{File.join(output_dir, 'gene.txt')}', 'wt');
      fprintf(fid, '%s\\n', gene{:});
      fclose(fid);

      fid = fopen('#{File.join(output_dir, 'mirna.txt')}', 'wt');
      fprintf(fid, '%s\\n', gene{:});
      fclose(fid);

      dlmwrite('#{File.join(output_dir, 'targets.txt')}', full(solucion));
    EOF

    pid, iin, iout, ierr = Open4.popen4(MATLAB)
    iin.write script
    iin.close
    iout.close
    ierr.close
    ignored, status = Process.waitpid2 pid
  end
end

if __FILE__ == $0
  endotelio = File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'data', 'examples', 'Endotelio')
  outputfile = File.join('/tmp/', 'milasso2')
  TaLasso.lasso(outputfile, endotelio, %w(tarbase))
end
