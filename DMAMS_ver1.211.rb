# DMAMS is a program written in Ruby ver. 2.4.
# The program detects genetic signature associated with selective sweep within the population.
#
# To run the script, the followings are required:
# 'fileutils'   # "gem install mechanize" for installation
# 'spreadsheet' # "gem install spreadsheet" for installation
# 'machanize'   # "gem install mechanize"
# 'daru'        # "gem install daru -v0.1.6" for installation
# 'statsample'  # "gem install statsample" for installation

# To run the DMAMS program, five arguments are required:
# 1) fasta_sequnce_file_name,
# 2) nwk_tree_fine_name,
# 3) "nuc" OR "pro" for clustering,
# 4) minimum_size_of_sequences_for_clustering, and
# 5) output_file_name.
#
# Example command: ruby DMAMS_ver1.211.rb example.fas example.nwk pro 0.05 test.csv
#
#
# Sequence file should be in single-line fasta format; the first sequence should be outgroup.
#
# Tree file should be in newick format; the tree must be rooted using the single outgroup.
#
# In the two input files, there must not be any symbols including "space" for sequence name except "underscore".
#
# Genetic signature can be either nucleotide or deduced amino acid. Please select the option either "nuc" or "pro".
# To detected genetic signature of deduced amino acid, sequence should be in-frame.
# (To detected genetic signature at single nucleotide, sequence can be untranslated region or out-frame coding region.)
#
# Minimum_size_of_strains_for_clustering is a proportion of strains of all sequences to determine subpopulation. The value should be smaller than 0.5
#
# Output file will be CSV.
# In the output, genetic signatures associated with significant selective sweep are marked '1' in the "detection" column.
# Using example files, there should be '1' for the row of position 10 amino acid L. The mutation was designated as a beneficial mutation to generate the sample data.
#
# Please contact Yuki FURUSE (furusey.kyoto@gmail.com) for any inquiry.


################################################################################
############ Please change file names and setting ##############################



fas_filename = ARGV[0]

tree_file_name = ARGV[1]

output_filename = ARGV[4]

if ARGV[2] == "nuc"
  nuc0_protein1 = 0 # 0 or 1, clustering based on...
elsif ARGV[2] == "pro"
  nuc0_protein1 = 1 # 0 or 1, clustering based on...
end  

cluster_size_temp = ARGV[3].to_f # minimum size for detection

############ You do not have to change anything below this line ################
################################################################################


print "DMAMS running"


prop_in_the_rest = 0.5 # smaller def 0.5
high_prop_initial = 0.9 # larger def 0.9

sig_p = 0.05 # default 0.05
sig_D = -1.5 # default -1.5


print "\n\nSequence file: ", fas_filename
print "\nTree file: ", tree_file_name
print "\nOutput file: ", output_filename, "\n"

if nuc0_protein1 == 0
  print "nucleotide clustering"
elsif nuc0_protein1 == 1
  print "protein clustering"
end

print "\nMinimum size of cluster to test: ", ARGV[3]


require 'fileutils'
#gem install mechanize
require 'daru' #gem install daru -v0.1.6
require 'statsample' #gem install statsample
mathE = Math::E

class Random
  include Math
    # Random number under Normal Distribution by Box-Muller
    def normal_rand(mu,sigma)
      a, b = self.rand(), self.rand()
      (sqrt(-2*log(rand()))*sin(2*PI*rand()) * sigma) + mu
    end

    # Random value under Poisson Distribution
    def poisson_rand(mu)
      lambda = Math.exp(-mu)
      k = 0
      p = 1.0
      while p >= lambda
        p *= self.rand()
        k += 1
      end
      return k - 1
    end
end

class Array
  def sum
    reduce(:+)
  end

  def mean
    sum.to_f / size
  end

  def var
    m = mean
    reduce(0) { |a,b| a + (b - m) ** 2 } / (size - 1)
  end

  def sd
    Math.sqrt(var)
  end
end



if nuc0_protein1 == 0
  atgc = ["A","T","G","C","-"]
else
  atgc = ["F", "L", "I", "M", "V", "S", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G", "X", "-", "B"]
end

codon_table_array =[["TTT", "F"],["TTC", "F"],["TTA", "L"],["TTG", "L"],
                    ["CTT", "L"],["CTC", "L"],["CTA", "L"],["CTG", "L"],
                    ["ATT", "I"],["ATC", "I"],["ATA", "I"],["ATG", "M"],
                    ["GTT", "V"],["GTC", "V"],["GTA", "V"],["GTG", "V"],
                    ["TCT", "S"],["TCC", "S"],["TCA", "S"],["TCG", "S"],
                    ["CCT", "P"],["CCC", "P"],["CCA", "P"],["CCG", "P"],
                    ["ACT", "T"],["ACC", "T"],["ACA", "T"],["ACG", "T"],
                    ["GCT", "A"],["GCC", "A"],["GCA", "A"],["GCG", "A"],
                    ["TAT", "Y"],["TAC", "Y"],["TAA", "B"],["TAG", "B"],
                    ["CAT", "H"],["CAC", "H"],["CAA", "Q"],["CAG", "Q"],
                    ["AAT", "N"],["AAC", "N"],["AAA", "K"],["AAG", "K"],
                    ["GAT", "D"],["GAC", "D"],["GAA", "E"],["GAG", "E"],
                    ["TGT", "C"],["TGC", "C"],["TGA", "B"],["TGG", "W"],
                    ["CGT", "R"],["CGC", "R"],["CGA", "R"],["CGG", "R"],
                    ["AGT", "S"],["AGC", "S"],["AGA", "R"],["AGG", "R"],
                    ["GGT", "G"],["GGC", "G"],["GGA", "G"],["GGG", "G"],["---", "-"]]



############# Tajima Definition ##############
$tajimaD = 0
def tajima(fas_seq)
  
  segregating_sites = 0
  difference = 0
  mean_pairwise_difference = 0
  a1 = 0
  a2 = 0
  b1 = 0
  b2 = 0
  c1 = 0
  c2 = 0
  e1 = 0
  e2 = 0
  $tajimaD = 0

  number_of_sequence = fas_seq.length 
  length_of_sequence = fas_seq[0].length
  
  i=0
  fas_seq[0].length.times do |i|
  
    check1 = fas_seq[0][i,1]

    fas_seq.each do |fas|
      check2 = fas[i,1]
      if check1 != check2
        segregating_sites += 1
        break
      end
    end
  end

  fas_seq.each do |fas1|
    fas_seq.each do |fas2|
      i=0
      length_of_sequence.times do |i|
        if fas1[i,1] != fas2[i,1]
          difference += 1
        end
      end
    end
  end

  mean_pairwise_difference = difference.to_f / (number_of_sequence * (number_of_sequence-1))

  (number_of_sequence-1).times do |j|
    a1 += 1.to_f/(j+1)
  end

  (number_of_sequence-1).times do |j|
    a2 += 1.to_f/((j+1)**2)
  end

  b1 = (number_of_sequence + 1).to_f / (3*(number_of_sequence - 1))

  b2 = (2 * ((number_of_sequence**2) + number_of_sequence + 3)).to_f / ((9 * number_of_sequence) * (number_of_sequence - 1))

  c1 = b1 - (1.to_f/a1)

  c2 = b2 - ((number_of_sequence + 2).to_f / (a1 * number_of_sequence)) + (a2.to_f / (a1**2))

  e1 = c1.to_f / a1

  e2 = c2.to_f / (a1**2 + a2)

  $tajimaD = (mean_pairwise_difference - (segregating_sites.to_f / a1)).to_f / (((e1 * segregating_sites) + (e2 * segregating_sites * (segregating_sites - 1)))**(1.to_f/2))
  
  if $tajimaD.finite?
  else
    $tajimaD = 999
  end
  
end # def end



############# make sequence array ##############
print "\n\nreading sequene\n"

fasta_array = []
temp_array = []
$number_of_sequence_original = 0
File.open(fas_filename) do |file|
  i=0
  file.each_line do |fas|
    temp_array << fas
    i += 1
    if i % 2 == 0
      fasta_array << temp_array
      temp_array = []
      $number_of_sequence_original += 1
    end
  end
end

cluster_size = ($number_of_sequence_original * cluster_size_temp).floor

if $number_of_sequence_original <= cluster_size * 2
  print "Number of sequences is too small OR Cluster size is too big"
  exit  ############# for single file/publication
end

if cluster_size < 10
  print "Number of sequences is too small OR Cluster size is too small"
  exit  ############# for single file/publication
end



high_prop = [1 - (cluster_size.to_f / $number_of_sequence_original), high_prop_initial].max.to_f



i=0 
fasta_array.each do |a_virus|
  a_virus.each do |fas|
    fas.slice!("\n")
    fas.slice!(">")
    i += 1
    if i % 2 == 0
      fas.upcase!
      fas.gsub!(/U/, 'T')
    end
  end
end

if nuc0_protein1 == 1
  fasta_array.each do |fas|
    fas[1] = fas[1].scan(/.{1,3}/)
  end
  fasta_array.each do |fas|
    temp_seq = ""
    fas[1].each do |codon_seq|
      good_protein = 0
      codon_table_array.each do |codon_ref|
        if codon_seq == codon_ref[0]
          temp_seq << codon_ref[1]
          good_protein = 1
          next
        end
      end
      if good_protein == 0
        temp_seq << "X"
      end
    end
    fas[1] = temp_seq
    temp_seq = ""
  end
end

root_strain = fasta_array[0][0]
$length_of_sequence_original = fasta_array[0][1].length



fasta_array_nuc = []
temp_array = []

File.open(fas_filename) do |file|
  i=0
  file.each_line do |fas|
    temp_array << fas
    i += 1
    if i % 2 == 0
      fasta_array_nuc << temp_array
      temp_array = []
    end
  end
end

i=0 
fasta_array_nuc.each do |a_virus|
  a_virus.each do |fas|
    fas.slice!("\n")
    fas.slice!(">")
    i += 1
    if i % 2 == 0
      fas.upcase!
      fas.gsub!(/U/, 'T')
    end
  end
end



########## Seg Site Array ################
fasta_array_temp = []
fasta_array.each do |fas|
  fasta_array_temp << fas[1]
end

seg_array = []
i = 0
$length_of_sequence_original.times do |i|
  check1 = fasta_array_temp[0][i,1]
  fasta_array_temp.each do |fas|
    check2 = fas[i,1]
    if check1 != check2
      seg_array  << i
      break
    end
  end
end



########## Read Tree File ###############
print "reading tree\n"


temp_str = ""
File.open(tree_file_name) do |file|
    i=0
    file.each_line do |nwk|
      if i == 0
        temp_str = nwk
        i += 1
      end
    end
end

temp_str.gsub!(/;/, "")
temp_str.slice!("\n")
temp_str.gsub!(/\(.*?/,"[")
temp_str.gsub!(/\).*?/,"]")

temp_array = []
tree_array = []

temp_array = temp_str.split("")

temp_str2 = ""


symbol1 = 1
symbol_brach_length = 0
temp_array.each do |text|
  if text == ":"
    symbol_brach_length = 1
    next
  end
  
  if symbol_brach_length == 1
    unless text == "[" || text == "]" || text == ","
      next
    end
  end
  
  symbol_brach_length = 0
      
  if text == "[" || text == "]" || text == ","
    symbol2 = 1
  else
    symbol2 = 0
  end
  if symbol1 == 1 && symbol2 == 0
   temp_str2 << '["'
  end

  if symbol1 == 0 && symbol2 == 1
   temp_str2 << '"]'
  end
  
  temp_str2 << text
  
  symbol1 = symbol2
end

 

tree_array = eval(temp_str2)
temp_array = []

if tree_array.length != 2
  print "Root?"
  exit
end



name_check = 0

tree_array.flatten.each do |strain_tree|
  fasta_array.each do |strain_fas|
    if strain_tree == strain_fas[0]
      name_check = 1
      break    
    end
  end
  
  if name_check == 0
    print "strain name not identical between files"
    exit
  end
  
  name_check = 0
end



######## Check each nuc distribution in entire tree #####
site_by_site_D = []

tajima_memory_array = []
strains_memory_array = []

print "checking genetic signature + calculating Tajima's D\n this may take a while\n\n"

print "progress...\n"

$length_of_sequence_original.times do |q|



  if q % 10 == 0
    print " ", ((q.to_f/$length_of_sequence_original)*100).round(1)
  end
  if q == $length_of_sequence_original - 1
    print " 100.0"
  end
  
  
  atgc.each do |nuc_acid|
    site_by_site_name = "#{q+1}#{nuc_acid}"
    
    lets_break = 0



########### Majority Minority Check ###############
lets_continue_majority = 0

on_number_in_all = 0
total_number_in_all = $number_of_sequence_original
i = 0

fasta_array.each do |fas|
  if fas[1][q,1] == nuc_acid
    on_number_in_all += 1
  end
end

if on_number_in_all.to_f <= total_number_in_all - (cluster_size * high_prop) && on_number_in_all >= cluster_size * high_prop
  lets_continue_majority = 1
end
    

    
##########################    
######## abortion ########
 if lets_continue_majority == 0
   next
 end
##########################



########## Divide Tree ###############

    tree_clusters = []
    
    tree_temp_array = []
    tree_temp_array = tree_array

    
    loop do  
      temp_array = []
      tree_temp_array_withNumber = []
      
      tree_temp_array.each do |sub_array|
        if sub_array.is_a?(String)
          temp_array << 1
        elsif sub_array.nil?
          temp_array << 0
        else
          temp_array << sub_array.flatten.length
        end
        temp_array << sub_array
        tree_temp_array_withNumber << temp_array
        temp_array = []
      end
    
    
      
      tree_temp_array_withNumber = tree_temp_array_withNumber.sort_by{ |aa, bb| aa }.reverse     # tree_temp_array_withNumber = [[number, big-subtree], [number, small-subtree]]
  
      tree_temp_array = []
      tree_temp_array_withNumber.each do |sub_array|
        tree_temp_array << sub_array[1]
      end

  
      if tree_temp_array_withNumber[0][0] < cluster_size || tree_temp_array_withNumber[0][0] < 4 || tree_temp_array_withNumber[0][0].nil?
        break
      end
      

      
      on_number_1 = 0
      
      fasta_array_temp2 = []
      fasta_array.each do |fas|
        if fas[1][q,1] == nuc_acid
          fasta_array_temp2 << fas[0]
        end
      end



      on_number_1 = (fasta_array_temp2 & tree_temp_array[0].flatten).length
 


      temp_array = []
      temp_array = tree_temp_array.drop(1)  # temp_array = tree of mid- and small-subtrees
      
      
      
      if on_number_1.to_f / tree_temp_array[0].flatten.length >= high_prop
        tree_clusters << tree_temp_array[0]
        tree_temp_array = temp_array
        
        
        if tree_temp_array.flatten.length < 2
          next
        end
      
      

      else
        temp_array = []  
        i = 0
        tree_temp_array.each do |sub_array|
          i += 1
          if i == 1 # divide the big-subtree
            temp_array << sub_array[0]
            temp_array << sub_array[1]
          else
            temp_array << sub_array # keep the small-subtree
          end
        end
        tree_temp_array = temp_array # now, tree contains 3 subtrees
        temp_array = []
      end
      


    end
    

    
############# Check prop in High and The rest ###################

    on_number_2 = 0
    the_rest_total = tree_temp_array.flatten.length
    
    fasta_array_temp2 = []
    fasta_array.each do |fas|
      if fas[1][q,1] == nuc_acid
        fasta_array_temp2 << fas[0]
      end
    end
    
    on_number_2 = (fasta_array_temp2 & tree_temp_array.flatten).length



##########################    
######## abortion ########
    if tree_clusters.flatten.nil? || tree_clusters.flatten.empty? || on_number_2.to_f / the_rest_total >= prop_in_the_rest || the_rest_total < cluster_size
      next
    end
##########################    



high_tajima = 0
low_tajima = 0
high_tajima_array = []
low_tajima_array = []



######## Tajima for High groups ########
    
   
    tree_clusters.each do |each_cluster|
      
      fasta_array_temp = []
      each_cluster.flatten.each do |strain|
        fasta_array_nuc.each do |fas|
          if fas[0] == strain
            fasta_array_temp << fas[1]
          end
        end
      end         
      
      $tajimaD = 0
      
      high_tajima = 0
      memory = 0
      strains_memory_array.each_with_index do |strains_set, array_index|
        if strains_set.flatten == each_cluster.flatten
          high_tajima = tajima_memory_array[array_index]
          memory = 1
          break
        end
      end
      
      if memory == 0
        tajima(fasta_array_temp)
        high_tajima = $tajimaD
        
        tajima_memory_array << high_tajima
        strains_memory_array << each_cluster.flatten
      end
      
      if each_cluster.flatten.length >= cluster_size
        high_tajima_array << high_tajima
      end


      
      each_cluster_temp = each_cluster
      
      temp_cluster = []
      loop do |yyy|
        if temp_cluster == each_cluster_temp.flatten(1)
          break
        end
        if yyy == 0
          temp_cluster = each_cluster_temp
        else
          temp_cluster = each_cluster_temp.flatten(1)
        end
        
        temp_cluster.each do |sub_cluster|
          if sub_cluster.kind_of?(Array)
            if sub_cluster.flatten.length >= 3
              fasta_array_temp = []
              sub_cluster.flatten.each do |strain|
                fasta_array_nuc.each do |fas|
                  if fas[0] == strain
                    fasta_array_temp << fas[1]
                  end
                end
              end



              $tajimaD = 0
      
              high_tajima_sub = 0
              memory = 0
              strains_memory_array.each_with_index do |strains_set, array_index|
                if strains_set.flatten == sub_cluster.flatten
                  high_tajima_sub = tajima_memory_array[array_index]
                  memory = 1
                  break
                end
              end
      
              if memory == 0
                tajima(fasta_array_temp)
                high_tajima_sub = $tajimaD
        
                tajima_memory_array << high_tajima_sub
                strains_memory_array << sub_cluster.flatten
              end

              if high_tajima_sub == 999
                next
              end


    
              if sub_cluster.flatten.length >= cluster_size
                high_tajima_array << high_tajima_sub
              end


              
            end
          end
        end
        each_cluster_temp = temp_cluster
      end
      




    
    ############# Tajima for Low group ###########
    
    



      
      fasta_array_temp = []
      tree_temp_array.flatten.each do |strain|
        if strain == root_strain
          next
        end
        fasta_array_nuc.each do |fas|
          if fas[0] == strain
            fasta_array_temp << fas[1]
          end
        end
      end
      
      
      $tajimaD = 0
      
      low_tajima = 0
      memory = 0
      strains_memory_array.each_with_index do |strains_set, array_index|
        if strains_set.flatten == tree_temp_array.flatten
          low_tajima = tajima_memory_array[array_index]
          memory = 1
        end
      end
      
      if memory == 0
        tajima(fasta_array_temp)
        low_tajima = $tajimaD
        
        tajima_memory_array << low_tajima
        strains_memory_array << tree_temp_array.flatten
      end
      
      if tree_temp_array.flatten.length >= cluster_size
        low_tajima_array << low_tajima
      end
      
      

############## Entire tree minus High #############


      
      each_cluster_temp = tree_array
      
      temp_cluster = []
      loop do |yyy|
        if temp_cluster == each_cluster_temp.flatten(1)
          break
        end
        
        temp_cluster = each_cluster_temp.flatten(1)
        
        sub_cluster_temp = []
        temp_cluster.each do |sub_cluster|
          if sub_cluster.kind_of?(Array)
            if sub_cluster.flatten.length >= 3
              fasta_array_temp_old = fasta_array_temp
              fasta_array_temp = []
              sub_cluster_temp = []
              sub_cluster.flatten.each do |strain|
                if strain == root_strain
                  next
                end
                
                temp_value = 0
                tree_clusters.flatten.each do |high_strain|
                  if strain == high_strain
                    temp_value = 1
                  end
                end
                if temp_value == 0
                  fasta_array_nuc.each do |fas|
                    if fas[0] == strain
                      fasta_array_temp << fas[1]
                      sub_cluster_temp << fas[0]
                    end
                  end
                end
                
                
                
              end
              
              if fasta_array_temp.length == 0 || fasta_array_temp_old == fasta_array_temp
                next
              end
              
              
              
              $tajimaD = 0
      
              low_tajima_sub = 0
              memory = 0
              strains_memory_array.each_with_index do |strains_set, array_index|
                if strains_set.flatten == sub_cluster_temp.flatten
                  low_tajima_sub = tajima_memory_array[array_index]
                  memory = 1
                end
               end
      
              if memory == 0
                tajima(fasta_array_temp)
                low_tajima_sub = $tajimaD
        
                tajima_memory_array << low_tajima_sub
                strains_memory_array << sub_cluster_temp.flatten
              end

              if low_tajima_sub == 999
                next
              end



              if fasta_array_temp.length >= cluster_size
                low_tajima_array << low_tajima_sub
              end

              

            end  # subcluster>3? end
          end  #subcluster array? end
        end  # temp_cluster each end
        each_cluster_temp = temp_cluster
      end  # loop_end



    
     
    
    temp_array = []
    temp_array << site_by_site_name
    temp_array << high_tajima
    temp_array << low_tajima
    temp_array << high_tajima_array
    
    
    
    temp_array << low_tajima_array
    
    
    temp_array_strains = []
    temp_array_strains << each_cluster.flatten
    
    temp_array << temp_array_strains
    
    
    site_by_site_D << temp_array
    temp_array = []
    temp_array_strains = []
    high_tajima_array = []
    low_tajima_array = []



    end # Tree Culster next ## for High groups
    


end # atgc end



end # whole length end



############## STAT ###########

print "\n\nstatistical test\n\n"
print "progress...\n"

site_by_site_diff = []

temp_array1 = []
temp_array2 = []
temp_array3 = []
temp_array4 = []


#site_by_site_D: 0-site, 1-highD, 2-restD, 3-highD_array, 4-restD_array, 5-high cluster strains (multi)



i = 0
site_by_site_D.each do |stat|
  i += 1
  if i % 10 == 0
    print " ", ((i.to_f/site_by_site_D.length)*100).round(1)
  end
  if i == site_by_site_D.length
    print " 100.0"
  end

  n1 = stat[3]
  n2 = stat[4]
  a1 = Daru::Vector[n1]
  b1 = Daru::Vector[n2]
  t1 = Statsample::Test::UMannWhitney.new(a1,b1)

  
  
  temp_array2 << i
  temp_array2 << stat[0]
  temp_array2 << stat[1]
  temp_array2 << stat[2]
  temp_array2 << (stat[1] - stat[2])
  temp_array2 << t1.probability_exact.to_f
  temp_array2 << stat[3]
  temp_array2 << stat[4]  
  temp_array2 << stat[5].flatten.length
  temp_array2 << stat[5].join(" / ")
  site_by_site_diff << temp_array2
  
  temp_array1 = []
  temp_array2 = []
end


#site_by_site_diff: 0-seqorder, 1-site, 2-highD, 3-restD, 4-diffD, 5-P,
                 # 6-highD array, 7-restD array, 8-high cluster number, 9-high cluter strains (join text)



  ###### BH correction #########
site_by_site_diff.each do |stat|
  if stat[2] < stat[3] && stat[2] < sig_D
    temp_array3 << stat[5]
    temp_array3 << stat[6]
    temp_array3 << stat[7]
    temp_array4 << temp_array3
    temp_array3 = []
  end
end

temp_array4.uniq!


temp_array4 = temp_array4.sort_by {|plist, highDarraylist, restDarraylist| plist}.reverse


i = 0
temp_array = []

length_of_list = temp_array4.length

value_bh = 0

temp_array4.each do |list|
  i += 1
  if list[0] < sig_p * ((length_of_list-i+1).to_f/length_of_list)
    value_bh = list[0]
    break
  end
end


#site_by_site_diff: 0-seqorder, 1-site, 2-highD, 3-restD, 4-diffD, 5-P,
                 # 6-highD array, 7-restD array, 8-high cluster number, 9-high cluter strains (join text)


  ###### Judgement ########

site_by_site_judge = []
temp_array = []

site_by_site_diff.each do |listitem|
  
  if listitem[2] == 999 || listitem[3] == 999
    temp_array << "."
  
   
  elsif listitem[5] <= value_bh && listitem[2] < listitem[3] && listitem[2] < sig_D
    temp_array << 1
  else
    temp_array << "."
  end
  
  temp_array << listitem[1][0...-1]
  temp_array << listitem[1][-1]
  temp_array << listitem[2].round(2)
  temp_array << listitem[3].round(2)
  temp_array << listitem[4].round(2)
  temp_array << listitem[5].round(4)
  temp_array << listitem[8]
  temp_array << listitem[9]
  
  site_by_site_judge << temp_array
  temp_array = []
end

#site_by_site_judge: 0-judgement, 1-site, 2-var, 3-highD, 4-restD, 5-diffD, 6-P, 7-high cluster number, 8-high cluter strains (join text)

temp_array = []
temp = ""
site_by_site_judge.each do |listitem_1|
  site_by_site_judge.each do |listitem_2|
    if listitem_1[8] == listitem_2[8]
      temp = listitem_2[1].to_s + listitem_2[2]
      temp_array << temp
      temp = ""
    end
  end
  listitem_1 << temp_array.join(" / ")
  temp_array = []
end

#site_by_site_judge: 0-judgement, 1-site, 2-var, 3-highD, 4-restD, 5-diffD, 6-P, 7-high cluster number, 8-high cluter strains (join text), 9-linked site

############# File output ###################################
temp_array = ["detection", "site", "var", "cluster_D", "rest_D", "diff", "Pvalue", "cluster size", "strains", "linked sites"]

site_by_site_judge.unshift(temp_array)


File.open("#{output_filename}", "a") do |file|
    site_by_site_judge.each do |list|
      i = 0
      list.each do |value|
        i +=1
          file.print(value)
        if i <= temp_array.length - 1
          file.print(",")
        end
      end
      file.print("\n")
    end
end



print "\n\ndone!\n"
