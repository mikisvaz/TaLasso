require 'open4'
module TaLasso
  MATLAB='matlab -nojvm -nodesktop -nodisplay'
  MATLAB_LIB= File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'Matlab')
  PUTATIVE_DIR = File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'data', 'targets')

  DEBUG=false

  def self.matlab_add_paths(lib_paths)
    script =<<-EOF
      addpath('#{MATLAB_LIB}');
#{
    lib_paths.collect {|path|
      File.join(MATLAB_LIB, path)
    }.collect {|path|
      "addpath('#{path}')"
    } * "\n"
}
    EOF
    script
  end

  def self.matlab_prepare_data(data_dir, putative)
    paths = putative.collect {|name| File.join(PUTATIVE_DIR, name) }

    script =<<-EOF
      [GE, ME, gene, mirna] = loadExpressionData('#{data_dir}');

      paths=[#{paths.collect {|path| "cellstr('" + path + "')" } * ', ' }];

      [GM, geneUnion, mirnaUnion] = loadPutativeData(paths);
      [GM, GEn, ME, gene, mirna]   = GMdataIntersect(GM, geneUnion, mirnaUnion, gene, mirna, GE, ME);

      GE = GEn - median(GEn, 2) * ones(1, size(GEn, 2));
    EOF
    script
  end
 
  def self.matlab_save(output_dir)
    script =<<-EOF
      fid = fopen('#{File.join(output_dir, 'gene.txt')}', 'wt');
      fprintf(fid, '%s\\n', gene{:});
      fclose(fid);

      fid = fopen('#{File.join(output_dir, 'mirna.txt')}', 'wt');
      fprintf(fid, '%s\\n', mirna{:});
      fclose(fid);

      dlmwrite('#{File.join(output_dir, 'targets.txt')}', full(solucion));
    EOF
    script
  end


  def self.matlab_run(script)
    pid, iin, iout, ierr = Open4.popen4(MATLAB)
    iin.write script
    iin.close
    if DEBUG
      ignored, status = Process.waitpid2 pid
      puts iout.read + ierr.read
      iout.close
      ierr.close
    else
      iout.close
      ierr.close
      ignored, status = Process.waitpid2 pid
    end
    status
  end

  def self.lasso(output_dir, data_dir, putative, lamb = 1/5)
    script  = matlab_add_paths(['miLasso', 'miLasso/l1_ls_matlab'])
    script += matlab_prepare_data(data_dir, putative)
    script +=<<-EOF
      solucion = l1_reg_model(GE, ME, GM, #{ lamb });
    EOF
    script += matlab_save(output_dir)

    FileUtils.mkdir_p output_dir unless File.exists? output_dir
    matlab_run(script)
  end

  def self.gen_mir(output_dir, data_dir, putative)
    script  = matlab_add_paths(['GenMiR'])
    script += matlab_prepare_data(data_dir, putative)
    script +=<<-EOF
      solucion = GenMiR(GE, ME, GM);
    EOF
    script += matlab_save(output_dir)

    FileUtils.mkdir_p output_dir unless File.exists? output_dir
    matlab_run(script)
  end

  def self.simple_gen_mir(output_dir, data_dir, putative)
    script  = matlab_add_paths(['SimpleGenMiR'])
    script += matlab_prepare_data(data_dir, putative)
    script +=<<-EOF
      solucion = GenMiRmodified(GE, ME, GM);
    EOF
    script += matlab_save(output_dir)

    FileUtils.mkdir_p output_dir unless File.exists? output_dir
    matlab_run(script)
  end


  def self.correlation(output_dir, data_dir, putative)
    script  = matlab_add_paths(['GenMiR'])
    script += matlab_prepare_data(data_dir, putative)
    script +=<<-EOF
      solucion = GEn./(std(GE')'*ones(1,size(GE,2)))*((ME-mean(ME,2)*ones(1,size(ME,2)))./(std(ME')'*ones(1,size(ME,2))))'/(size(GE,2)-1).*GM;
    EOF
    script += matlab_save(output_dir)

    FileUtils.mkdir_p output_dir unless File.exists? output_dir
    matlab_run(script)
  end


end

if __FILE__ == $0
  endotelio = File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'data', 'examples', 'Endotelio')
  outputfile = File.join('/tmp/', 'milasso2')
  TaLasso.lasso(outputfile, endotelio, %w(tarbase))
end
