module MatrixFormat
  def self.matrix2list(file_matrix, file_row_names, file_col_names, sep=',')
    row_names  = File.open(file_row_names).read.split("\n")
    col_names = File.open(file_col_names).read.split("\n")

    list = {}
    File.open(file_matrix).read.split("\n").each_with_index do |line, i|
      row_name = row_names[i]
      values = line.chomp.split(sep)
      col_names.each_with_index do |col_name, i|
        value = values[i].to_f
        next if value == 0
        list[row_name] ||= {}
        list[row_name][col_name] =  value
      end
    end

    list
  end

  def self.read_list(file_list, sep="\t")
    list = {}
    File.open(file_list).each do |line|
      row, col, value = line.chomp.split("\t")
      list[row] = {}
      list[row][col] = value || 1
    end
    list 
  end

  def self.list2matrix(list, file_matrix, file_row_names, file_col_names, sep =",")
    row_names = list.keys
    col_names = list.collect {|key, values| values.keys }.flatten.uniq


    File.open(file_row_names, 'w') do |f| f.puts row_names * "\n" end
    File.open(file_col_names, 'w') do |f| f.puts col_names * "\n" end
    File.open(file_matrix, 'w') do |f| 
      f.puts row_names.collect {|row_name|
        info = list[row_names] || {}
        col_names.collect { |col_name|
          info[col_name] || 0
        } * sep
      } * "\n"
    end
  end
end
