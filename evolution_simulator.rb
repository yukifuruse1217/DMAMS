print "evolution.rb running\n"

require 'csv'
require 'fileutils'

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



$amino_acid = "empty"
$amino_acid_child = "empty"
$amino_acid_parent = "empty"
def translation(a,b,c)
	if [a,b,c] == [1,1,1] || [a,b,c] == [1,1,3]
		$amino_acid = "Phe"
	elsif [a,b,c] == [1,1,0] || [a,b,c] == [1,1,2] || [a,b,c] == [3,1,0] || [a,b,c] == [3,1,1] || [a,b,c] == [3,1,2] || [a,b,c] == [3,1,3]
		$amino_acid = "Leu"
	elsif [a,b,c] == [0,1,0] || [a,b,c] == [0,1,1] || [a,b,c] == [0,1,3]
		$amino_acid = "Ile"
	elsif [a,b,c] == [0,1,2]
		$amino_acid = "Met"
	elsif [a,b,c] == [2,1,0] || [a,b,c] == [2,1,1] || [a,b,c] == [2,1,2] || [a,b,c] == [2,1,3]
		$amino_acid = "Val"
	elsif [a,b,c] == [1,3,0] || [a,b,c] == [1,3,1] || [a,b,c] == [1,3,2] || [a,b,c] == [1,3,3] || [a,b,c] == [0,2,1] || [a,b,c] == [0,2,3]
		$amino_acid = "Ser"
	elsif [a,b,c] == [3,3,0] || [a,b,c] == [3,3,1] || [a,b,c] == [3,3,2] || [a,b,c] == [3,3,3]
		$amino_acid = "Pro"
	elsif [a,b,c] == [0,3,0] || [a,b,c] == [0,3,1] || [a,b,c] == [0,3,2] || [a,b,c] == [0,3,3]
		$amino_acid = "Thr"
	elsif [a,b,c] == [2,3,0] || [a,b,c] == [2,3,1] || [a,b,c] == [2,3,2] || [a,b,c] == [2,3,3]
		$amino_acid = "Ala"
	elsif [a,b,c] == [1,0,1] || [a,b,c] == [1,0,3]
		$amino_acid = "Tyr"
	elsif [a,b,c] == [3,0,1] || [a,b,c] == [3,0,3]
		$amino_acid = "His"
	elsif [a,b,c] == [3,0,0] || [a,b,c] == [3,0,2]
		$amino_acid = "Gln"
	elsif [a,b,c] == [0,0,1] || [a,b,c] == [0,0,3]
		$amino_acid = "Asn"
	elsif [a,b,c] == [0,0,0] || [a,b,c] == [0,0,2]
		$amino_acid = "Lys"
	elsif [a,b,c] == [2,0,1] || [a,b,c] == [2,0,3]
		$amino_acid = "Asp"
	elsif [a,b,c] == [2,0,0] || [a,b,c] == [2,0,2]
		$amino_acid = "Glu"
	elsif [a,b,c] == [1,2,1] || [a,b,c] == [1,2,3]
		$amino_acid = "Cys"
	elsif [a,b,c] == [1,2,2]
		$amino_acid = "Trp"
	elsif [a,b,c] == [3,2,0] || [a,b,c] == [3,2,1] || [a,b,c] == [3,2,2] || [a,b,c] == [3,2,3] || [a,b,c] == [0,2,0] || [a,b,c] == [0,2,2]
		$amino_acid = "Arg"
	elsif [a,b,c] == [2,2,0] || [a,b,c] == [2,2,1] || [a,b,c] == [2,2,2] || [a,b,c] == [2,2,3]
		$amino_acid = "Gly"
	elsif [a,b,c] == [1,0,0] || [a,b,c] == [1,0,2] || [a,b,c] == [1,2,0]
		$amino_acid = "stop"
	end
end















######## SETTING ########
output_file = "evo19_neut_popLinDec"

sweep0_diversify1_conv2_neut3 = 3

simulation_number = 100



default_pickup = 1000 # def 1000
if sweep0_diversify1_conv2_neut3 == 0
  good_degree = 2000 # sweep def 2000
elsif sweep0_diversify1_conv2_neut3 == 1
  good_degree = 4000 # div def 4000
  bad_degree = 100 # div def 100, to decrease synonymous mutation
elsif sweep0_diversify1_conv2_neut3 == 2
  good_degree = 4000 # conv def 4000
  bad_degree = 100   # conv def 100, to decrease good-to-bad mutation
end



population_size = 1000 # def 200 (100/400,   10 for increase, 1000 for decrease)

pop_no0_plus1_multi2 = 1 # 0, 1 or 2
pop_change = -0.8  # def 0.0, can be negative to decrease (0.39 for plus, 0.0035 for multi, -0.8 for minus, -0.0016 for div)

no_mutation_for_the_first = 200 # def 1000 no mutation only at good position;   when pop400, 2000;   when neut && pop=400/Lin/Exp/Decline end-1000
#no_good_effect_for_the_first = 10 # def 10


end_generation = 1200 # def 2000 for neut/div;  when pop=Lin/Exp, 1100;  when pop=400, 2100;   for Pop decline neut, 1200
#end_Ala_pop = 0.5 #def 0.5
end_good_population = 0.3 # def 0.3 for conv/sweep


codon_length = 1000 # def 1000
sample_number = 200 # def 200
sample_number_pro = 50 # def 50
mutation_per_site_generation = 0.0002 # def 0.0002

if sweep0_diversify1_conv2_neut3 == 2
  mutation_per_site_generation_good = 0.002 # def  0.002 for conv
else
  mutation_per_site_generation_good = 0.002 # def 0.002 for neut/sweep/div
end

# A=0, T/U=1, G=2, C=3:
initial_seq = [2,3,1] # def 231
# Ala=GCU
# AAA=Asn, UUU=Phe, GGG=Gly, CCC=Pro
#good_codonset = [[3,1,0],[3,1,1],[3,1,2],[3,1,3]]
# Leu: CUX=31X
# Thr: ACX=03X
good_aminoacid = "Leu" # def "Leu"
good_position = 9 # def 9 = position 10
#good_term = (0..100).to_a    ## CHECK ON/OFF!
sample_freq = 100 # dif 100
######## SETTING ########


######## simulation repeat########
j = 0
simulation_number.times do
  j += 1
  
  
#  j = 36
  
  
  print "\n", "simulation #", j, "\n"
######## simulation repeat########







$t_generation = 0
nonsense_codonset = [[1,0,0],[1,0,2],[1,2,0]]
#UAA=100, UAG=102, UGA=120
nonsense_check = 0
a_virus = []
viruses = []
viruses_now = []
viruses_sampled = []
new_codon = []
new_nucleotide = 0
#progeny = 0
progeny_array = []
population_size.times do
  progeny_array << 0
end
progeny_array_new = []

############# START ################


a_virus=[]


codon_length.times do
  a_virus << initial_seq
end

population_size.times do
  viruses << a_virus
end


############# next generation begin ############                 
q = 0
loop do
  q += 1

  
  viruses_now=[]
  viruses_for_progeny = []
  progeny_array_new = []
  temp_array = []
  

  viruses.length.times do |i|
    if q == 1
      default_pickup.times do
        temp_array << viruses[i]
      end
    
    else
      
      
      if progeny_array[i] == 0
        default_pickup.times do
          temp_array << viruses[i]
        end
      
      elsif progeny_array[i] == 1
        good_degree.times do
          temp_array << viruses[i]
        end
      elsif progeny_array[i] == 2
        bad_degree.times do
          temp_array << viruses[i]
        end
      
      else
        print "progeny_array?"
      end
    
    end
  end
  
  
  if pop_no0_plus1_multi2 == 0 
    population_size.times do
    viruses_for_progeny << temp_array.sample(1)
  end
  end

  if pop_no0_plus1_multi2 == 1 
  [50,((population_size.to_f+((q-1)*pop_change).floor).to_f)].max.round.times do
    viruses_for_progeny << temp_array.sample(1)
  end
  end

  if pop_no0_plus1_multi2 == 2
  [50,(((population_size.to_f*((1+pop_change)**(q-1))).floor).to_f)].max.round.times do
    viruses_for_progeny << temp_array.sample(1)
  end
  end
  


#  print " ", viruses_for_progeny.length, " "  
  
  
  viruses_for_progeny = viruses_for_progeny.flatten(1)
  
  
  
    
  viruses_for_progeny.each do |strain|
    

############# generate mutation ############                 
          nonsense_check = 1
          
          until nonsense_check == 0 do
            a_virus=[]
            strain.each_with_index do |triplet_codon, position_codon|
              new_codon=[]
              
              
              if q <= no_mutation_for_the_first && position_codon == good_position # No mutation at position 10 for the fist gen
                triplet_codon.each do |single_nucleotide|
                  new_nucleotide = single_nucleotide
                  new_codon << new_nucleotide
                end
              
              else
                triplet_codon.each do |single_nucleotide|
                  if position_codon == good_position
                    mutation_rate = mutation_per_site_generation_good
                  else
                    mutation_rate = mutation_per_site_generation
                  end
                                    
                  if Random.rand > mutation_rate
                    new_nucleotide = single_nucleotide
                  else
                    loop do
                      new_nucleotide = Random.rand(4)
                      if new_nucleotide != single_nucleotide
                        break
                      end
                    end
                  end
                  new_codon << new_nucleotide
                end
              
              end
              a_virus << new_codon
            end

############# go back to untill, if nonsense ############                 
            nonsense_check = 0
            a_virus.each do |check_codon|
              nonsense_codonset.each do |nonsense_codon|
                if check_codon == nonsense_codon
                  nonsense_check = 1
                end
              end
            end
            
#           if q <= no_mutation_for_the_first # No mutation at position 10 for the fist gen
#             $amino_acid = "empty"
#             translation(*a_virus[good_position])
#               if $amino_acid != "Ala"
#                 nonsense_check = 1
#               end
#           end

#            if q <= 10 # No GOOD mutation at position 10 for the fist 10 gen
#              $amino_acid = "empty"
#              translation(*a_virus[good_position])
#              if sweep0_diversify1_conv2_neut3 == 1
#                if $amino_acid != "Ala"
#                  nonsense_check = 1
#                end
#              elsif sweep0_diversify1_conv2_neut3 == 0 || sweep0_diversify1_conv2_neut3 == 2
#                if $amino_acid == good_aminoacid
#                  nonsense_check = 1
#                end
#              end
#            end
                        
            
          end



############# good check ############                 
            
if sweep0_diversify1_conv2_neut3 == 0
  $amino_acid = "empty"
  $amino_acid_child = "empty"
    
  translation(*a_virus[good_position])
  $amino_acid_child = $amino_acid
  $amino_acid = "empty"
            
  if $amino_acid_child == good_aminoacid
    progeny_array_new << 1
  else
    progeny_array_new << 0
  end
  
  
  
elsif sweep0_diversify1_conv2_neut3 == 1
  $amino_acid = "empty"
  $amino_acid_child = "empty"
  $amino_acid_parent = "empty"
  
  translation(*a_virus[good_position])
  $amino_acid_child = $amino_acid
  $amino_acid = "empty"
            
  translation(*strain[good_position])
  $amino_acid_parent = $amino_acid
  $amino_acid = "empty"
  
  
  if $amino_acid_child == $amino_acid_parent
    if a_virus[good_position] == strain[good_position]
      progeny_array_new << 0
    else
      progeny_array_new << 2
    end
  else
    progeny_array_new << 1
  end

elsif sweep0_diversify1_conv2_neut3 == 2
  $amino_acid = "empty"
  $amino_acid_child = "empty"
  $amino_acid_parent = "empty"
  
  translation(*a_virus[good_position])
  $amino_acid_child = $amino_acid
  $amino_acid = "empty"
            
  translation(*strain[good_position])
  $amino_acid_parent = $amino_acid
  $amino_acid = "empty"
            
  if $amino_acid_child != $amino_acid_parent && $amino_acid_child == good_aminoacid 
    progeny_array_new << 1
  elsif $amino_acid_child != $amino_acid_parent && $amino_acid_parent == good_aminoacid
    progeny_array_new << 2
  else
    progeny_array_new << 0
  end

elsif sweep0_diversify1_conv2_neut3 == 3
  progeny_array_new << 0
end
            











############# determine set of this generation ############                         
        
          viruses_now << a_virus
 
  end

  viruses=[]
  viruses = viruses_now
  progeny_array = []
  progeny_array = progeny_array_new
  $t_generation += 1

########################## display for check
print("Gen ", $t_generation, "   ")
#print(viruses.length, "   ")
#viruses.each do  |strain|
#  print strain, "\n"
#end
########################## display for check






########################## file output
if $t_generation % sample_freq == 0
File.open("#{output_file}-seq#{j}.fas","a") do |file|
    
      i = 0
      if $t_generation == sample_freq
        file.print(">0_Ala", "\n")
        codon_length.times do
          initial_seq.each do |seq_num|
            if seq_num == 0
              seq_nuc = "A"
            elsif seq_num == 1
              seq_nuc = "T"
            elsif seq_num == 2
              seq_nuc = "G"
            elsif seq_num == 3
            seq_nuc = "C"
            end
          file.print(seq_nuc)
          end
        end
        file.print("\n")
      end
  

             ######## ALL or random samples #########  
     viruses_sampled = []
      if viruses.length < sample_number
        viruses_sampled = viruses
      else
        viruses_sampled = viruses.sample(sample_number)
      end
      viruses_sampled.each do |strain|
        i += 1
        file.print(">", $t_generation, "_", i)
        
        
        $amino_acid = "empty"
        translation(*strain[good_position])
        if $amino_acid == good_aminoacid
          file.print("_", $amino_acid, "_Good")
        else
          file.print("_", $amino_acid)
        end
        
        

        
        file.print("\n")
        
        
        
        strain.flatten.each do |seq_num|
          if seq_num == 0
            seq_nuc = "A"
          elsif seq_num == 1
            seq_nuc = "T"
          elsif seq_num == 2
            seq_nuc = "G"
          elsif seq_num == 3
          seq_nuc = "C"
          end
          file.print(seq_nuc)
        end
        file.print("\n")
      end
             ######## samples #########      
end





File.open("#{output_file}-pop.csv","a") do |file|
    i = 0
    ii = 0
    
    if j == 1 && $t_generation == sample_freq
      file.print("sim,gen,total_pop,Ala_pop,ALA_prop,good_pop,good_prop", "\n")
    end    
    
    
    viruses.each do |strain|
      translation(*strain[good_position])
      if $amino_acid == "Ala"
        i += 1
      end
      if $amino_acid == good_aminoacid
        ii += 1
      end
    end
    file.print(j, ",", $t_generation, ",", viruses.length, ",", i, ",", ((i.to_f/viruses.length)*100).round(2), ",", ii, ",", ((ii.to_f/viruses.length)*100).round(2), "\n")
end



end # record frequency
########################## file output  





############# stop generation criteria ############ 

break_check = 0

if sweep0_diversify1_conv2_neut3 == 0 ||  sweep0_diversify1_conv2_neut3 == 2
  temp_var = 0
  viruses.each do |strain|
    $amino_acid = "empty"
    translation(*strain[good_position])
    if $amino_acid == good_aminoacid
      temp_var += 1
    end
  end
  
  if temp_var > viruses.length.to_f * end_good_population
    break_check = 1
  end


#if sweep0_diversify1_conv2_neut3 == 0 || sweep0_diversify1_conv2_neut3 == 1 || sweep0_diversify1_conv2_neut3 == 2
#  temp_var = 0
#  viruses.each do |strain|
#    $amino_acid = "empty"
#    translation(*strain[good_position])
#    if $amino_acid == "Ala"
#      temp_var += 1
#    end
#  end
  
#  if temp_var < viruses.length * end_Ala_pop
#    break_check = 1
#  end





  
elsif sweep0_diversify1_conv2_neut3 == 1 || sweep0_diversify1_conv2_neut3 == 3
  if end_generation <= $t_generation
    break_check = 1
  end
end


#if end_generation <= $t_generation
#  break
#end






################ Last Record #################
if break_check == 1

print "\n\nfasta_edit running\n"



viruses_sampled = viruses.sample(sample_number)




if $t_generation % sample_freq != 0
File.open("#{output_file}-pop.csv","a") do |file|
    i = 0
    ii = 0
    
    
    viruses.each do |strain|
      translation(*strain[good_position])
      if $amino_acid == "Ala"
        i += 1
      end
      if $amino_acid == good_aminoacid
        ii += 1
      end
    end
    file.print(j, ",", $t_generation, ",", viruses.length, ",", i, ",", ((i.to_f/viruses.length)*100).round(2), ",", ii, ",", ((ii.to_f/viruses.length)*100).round(2), "\n")
end
end






#if $t_generation % sample_freq == 0
File.open("#{output_file}-seq#{j}-0andLast.fas","a") do |file|
    
      i = 0
#      if $t_generation == 1
        file.print(">0_Ala", "\n")
        codon_length.times do
          initial_seq.each do |seq_num|
            if seq_num == 0
              seq_nuc = "A"
            elsif seq_num == 1
              seq_nuc = "T"
            elsif seq_num == 2
              seq_nuc = "G"
            elsif seq_num == 3
            seq_nuc = "C"
            end
          file.print(seq_nuc)
          end
        end
        file.print("\n")
#      end
  



             ######## random samples #########  
      viruses_sampled.each do |strain|
        i += 1
        file.print(">", $t_generation, "_", i)
        
        
        $amino_acid = "empty"
        translation(*strain[good_position])
        if $amino_acid == good_aminoacid
          file.print("_", $amino_acid, "_Good")
        else
          file.print("_", $amino_acid)
        end
        


        
        file.print("\n")
        
        
        
        strain.flatten.each do |seq_num|
          if seq_num == 0
            seq_nuc = "A"
          elsif seq_num == 1
            seq_nuc = "T"
          elsif seq_num == 2
            seq_nuc = "G"
          elsif seq_num == 3
          seq_nuc = "C"
          end
          file.print(seq_nuc)
        end
        file.print("\n")
      end
             ######## random samples #########      
end











viruses_sampled_protein = viruses_sampled.sample(sample_number_pro)





#if $t_generation % sample_freq == 0
File.open("#{output_file}-seq#{j}-0andLast_sample.fas","a") do |file|
    
      i = 0
#      if $t_generation == 1
        file.print(">0_Ala", "\n")
        codon_length.times do
          initial_seq.each do |seq_num|
            if seq_num == 0
              seq_nuc = "A"
            elsif seq_num == 1
              seq_nuc = "T"
            elsif seq_num == 2
              seq_nuc = "G"
            elsif seq_num == 3
            seq_nuc = "C"
            end
          file.print(seq_nuc)
          end
        end
        file.print("\n")
#      end
  



             ######## random samples #########  
      viruses_sampled_protein.each do |strain|
        i += 1
        file.print(">", $t_generation, "_", i)
        
        
        $amino_acid = "empty"
        translation(*strain[good_position])
        if $amino_acid == good_aminoacid
          file.print("_", $amino_acid, "_Good")
        else
          file.print("_", $amino_acid)
        end
        


        
        file.print("\n")
        
        
        
        strain.flatten.each do |seq_num|
          if seq_num == 0
            seq_nuc = "A"
          elsif seq_num == 1
            seq_nuc = "T"
          elsif seq_num == 2
            seq_nuc = "G"
          elsif seq_num == 3
          seq_nuc = "C"
          end
          file.print(seq_nuc)
        end
        file.print("\n")
      end
             ######## random samples #########      
end























#end







############# FasToMeg #############

fas_filename = "#{output_file}-seq#{j}-0andLast.fas"
#     #{ii+1}



############# make sequence array ##############
name_array = []
fasta_array = []

temp_array = []
number_of_sequence_original = 0

File.open(fas_filename) do |file|
  i=0
  file.each_line do |fas|
    temp_array = fas
    temp_array.slice!("\n")
    temp_array.slice!(">")
    i += 1
    if i % 2 == 0
      fasta_array << temp_array
      temp_array = []
      number_of_sequence_original += 1
    else
      name_array << temp_array
      temp_array = []
    end
  end
end


if fasta_array.length != name_array.length
  print "\n fasta file line numbers error"
end

nchar = 0
fasta_array.each do |fas|
  if fas.length > nchar
    nchar = fas.length
  end
end


############# make MEGA file ##############
File.open("#{fas_filename}.meg","a") do |file|
    file.print(
      "#mega", "\n",
      "!Title ;", "\n",
      "!Format DataType=DNA indel=- CodeTable=Standard;", "\n",
      "\n",     
      "!Domain=Data property=Coding CodonStart=1;", "\n")

    fasta_array.length.times do |i|
      file.print("#", name_array[i], "\n", fasta_array[i], "\n")
    end
    
end
               
               

############### FatToMega END ###############




































############# FasToMeg #############

fas_filename = "#{output_file}-seq#{j}-0andLast_sample.fas"
#     #{ii+1}



############# make sequence array ##############
name_array = []
fasta_array = []

temp_array = []
number_of_sequence_original = 0

File.open(fas_filename) do |file|
  i=0
  file.each_line do |fas|
    temp_array = fas
    temp_array.slice!("\n")
    temp_array.slice!(">")
    i += 1
    if i % 2 == 0
      fasta_array << temp_array
      temp_array = []
      number_of_sequence_original += 1
    else
      name_array << temp_array
      temp_array = []
    end
  end
end


if fasta_array.length != name_array.length
  print "\n fasta file line numbers error"
end

nchar = 0
fasta_array.each do |fas|
  if fas.length > nchar
    nchar = fas.length
  end
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
                    ["GGT", "G"],["GGC", "G"],["GGA", "G"],["GGG", "G"]]




fasta_array = []
temp_array = []
File.open("#{output_file}-seq#{j}-0andLast_sample.fas") do |file|
  i=0
  file.each_line do |fas|
    temp_array << fas
    i += 1
    if i % 2 == 0
      fasta_array << temp_array
      temp_array = []
    end
  end
end

i=0 
fasta_array.each do |a_virus|
  a_virus.each do |fas|
    fas.slice!("\n")
  end
end


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

fasta_array.flatten.each do |fas|
  File.open("#{output_file}-seq#{j}-0andLast_sample-Protein.fas","a") do |filefile|
    filefile.print(fas)
    filefile.print("\n")
  end
end

































fas_filename = "#{output_file}-seq#{j}-0andLast_sample-Protein.fas"



############# make sequence array ##############
name_array = []
fasta_array = []

temp_array = []
number_of_sequence_original = 0

File.open(fas_filename) do |file|
  i=0
  file.each_line do |fas|
    temp_array = fas
    temp_array.slice!("\n")
    temp_array.slice!(">")
    i += 1
    if i % 2 == 0
      fasta_array << temp_array
      temp_array = []
      number_of_sequence_original += 1
    else
      name_array << temp_array
      temp_array = []
    end
  end
end


if fasta_array.length != name_array.length
  print "\n fasta file line numbers error"
end

nchar = 0
fasta_array.each do |fas|
  if fas.length > nchar
    nchar = fas.length
  end
end


############# make nex file ##############
File.open("#{fas_filename}.nxs","a") do |file|
    file.print(
      "#NEXUS", "\n",
      "begin data;", "\n",
    	"\t", "dimensions ntax=", fasta_array.length, " nchar=", nchar, ";", "\n",
    	"\t", "format datatype=protein missing=? gap=-;", "\t", "\n",
      "matrix", "\n")

    fasta_array.length.times do |i|
      file.print(name_array[i], "     ", fasta_array[i], "\n")
    end
    
    file.print(
    ";", "\n",
    "end;")

end
               
               


############### FasToNex  END  #####################











































  break #next simulation
end ### Record on off

######LAST RECORD ########################################################


























































############# next generation repeat ############ 
end
############# next generation repeat ############ 





######## simulation repeat########
end
######## simulation repeat########




























































print "\ndone!\n"
