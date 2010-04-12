$LOAD_PATH.unshift(File.join(File.dirname(__FILE__), '..', 'lib'))
require 'matrix_format'
require 'inline'
require 'gruff'
require 'MARQ/score'

module TaLasso
  module TaLasso::Validation

  class << self
      inline do |builder|
        builder.c_raw <<-EOC
  
    /**
    * Compute log(k!)
    * @param k The value k.
      * @return The result.
      */
      double lFactorial(double k)
    {
      double r = 0;
      int i;
      for(i=2 ; i<=(int)k ; i++)
        {
          r = r + (double)(log((double)i));
        }
        return r;
    }
  
  
  
    /**
      * Compute the log(binom(n,k))
      * @param n The number of possible items.
      * @param k The number of selected items.
      * @return The result.
      */
      double lBinom(double n, double k)
    {
      long i;
      double r = 0;
      if(n > n-k){
        k = n-k;
      }
  
      for(i = (long)n ; i> (n-k) ; i--)
        {
          r = r + log((double)i);
        }
        r = r - lFactorial(k);
        return r;
    }
    EOC
  
    builder.c <<-EOC
    /**
      * Compute the Hypergeometric accumulated value.
      * @param total => total size
      * @param support => total support
      * @param list => selected list size,
      * @param found => support
      * @return The result
      */
      double hypergeometric(double total, double support, double list, double found)
    {
      double other = total - support;
  
      double top = list;
      if(support < list){
        top = support;
      }
  
      double log_n_choose_k = lBinom(total,list);
  
      double lfoo = lBinom(support,top) + lBinom(other, list-top);
      double sum = 0;
      int i;
      for (i = (int)top; i >= found; i-- )
        {
          sum = sum + exp(lfoo - log_n_choose_k);
          if ( i > found)
            {
              lfoo = lfoo + log(i / (support - i+1)) + log( (other - list + i) / (list-i+1) );
            }
        }
        return sum;
    }
    EOC
      end
    end

    def self.draw_pvalues(pvalues, size, filename)
      require 'gnuplot'
      sizes, values = pvalues.collect.sort_by{|p| p[0]}.transpose.values_at(0, 1)
      Gnuplot.open do |gp|
        Gnuplot::Plot.new( gp ) do |plot|
          plot.title  "P-value distribution"
          plot.xrange "[0:#{size}]"
          plot.yrange "[0:1]"

          plot.ylabel "p-value"
          plot.xlabel "Size"

          plot.terminal "png"
          plot.output filename

          plot.data << Gnuplot::DataSet.new( [sizes, values] ) do |ds|
            ds.with = "lines"
            ds.linewidth = 1
          end
        end
      end
    end

    def self.positions(results_path, db_path)
      list = MatrixFormat.matrix2list(File.join(results_path, 'targets.txt'), File.join(results_path, 'gene.txt'), File.join(results_path, 'mirna.txt'))
      gs   = MatrixFormat.matrix2list(File.join(db_path, 'putative.txt'), File.join(db_path, 'gene.txt'), File.join(db_path, 'mirna.txt'))

      result_pairs = list.collect{|gene, targets|
        targets.collect{|mirna, value|
          {:gene => gene, :mirna => mirna, :value => value.to_f}
        }
      }.flatten.sort_by{|info| info[:value]}.collect{|info| [info[:gene], info[:mirna]] * ":"}

      gs_pairs = gs.collect{|gene, targets|
        targets.collect{|mirna, value|
          {:gene => gene, :mirna => mirna, :value => value}
        }
      }.flatten.collect{|info| [info[:gene], info[:mirna]] * ":"}

      positions = gs_pairs.collect{|pair| result_pairs.index(pair)}.compact.sort

      {:results => result_pairs.length, :positions => positions}
    end

    def self.scale(positions, results, size = 10000.0)
      ratio = results.to_f / size

      positions.collect{|pos| (pos / ratio).to_i }
    end

    def self.validate(positions, results, filename = nil, steps = 1000)
      sizes = (1..steps).collect{|i| (i * results / steps).to_i }
      goldstandard = positions.length

      pvalues = {}
      sizes.collect do |total_targets|
        next if total_targets > results
        pvalues[total_targets] = hypergeometric(results, goldstandard, total_targets, positions.count_smaller(total_targets))
      end

      draw_pvalues(pvalues, results, filename) if filename
      pvalues
    end

    def self.hit_score(positions, results, filename = nil)
      max         = 10000

      if results > max
        positions   = scale(positions, results, max)
        results     = max
      end

      score       = Score.fast_score_scale(positions, results, 0)

      times       = 5000
      null_scores = Score.null_scores(results, 0, times).sort

      pvalue      = (times - null_scores.count_smaller(score.abs)).to_f / times 

      Score.draw_hits(positions, results, filename, :size => 2000, :bg_color => :green ) if filename
      {:score => score, :pvalue => pvalue}
    end
  end
end

if __FILE__ == $0
  PUTATIVE_DIR = File.join(File.dirname(File.dirname(File.expand_path(__FILE__))), 'data', 'targets')
  info = TaLasso::Validation.positions('/tmp/milasso2',  File.join(PUTATIVE_DIR, 'tarbase'))
  p info
  p TaLasso::Validation.validate(info[:positions], info[:results], '/tmp/pvalues.png').values.sort
  p TaLasso::Validation.hit_score(info[:positions], info[:results], '/tmp/hits.png')
end
