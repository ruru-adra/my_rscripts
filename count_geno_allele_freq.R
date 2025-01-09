#input data must in long and filter alel=0/no alel detected


#Calculate genotype frequencies
genotype_freq<- fltr_bk_boer_alel %>%
  group_by(snpID, alel, Breed) %>%  # Use 'alel' as is
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(snpID, Breed) %>%
  mutate(Frequency = round(Count / sum(Count), 3)) %>%
  ungroup() %>%
  # Calculate frequencies for the total population
  group_by(snpID, alel) %>%
  mutate(Total_Population_Frequency = sum(Count) / sum(Count[!duplicated(Breed)])) %>%
  select(snpID, alel, Breed, Frequency, Total_Population_Frequency, Count)

# Reformat the data to match the desired format
genotype_freq_wide <- genotype_freq %>%
  pivot_wider(names_from=Breed, values_from=c(Frequency, Count)) %>%
  arrange(snpID, alel)

#Calculate allele freq
#Split the Genotype into two Alleles
allele_freq <- geno_data %>%
  mutate(Allele1 = substr(alel, 1, 1), Allele2 = substr(alel, 2, 2)) %>%
  # Gather to long format to treat both alleles equally
  gather(key = "Allele_Position", value = "Allele", Allele1, Allele2) %>%
  group_by(snpID, Allele) %>%
  # Calculate allele frequencies for each breed
  summarise(
    Frequency_BK = sum(Frequency_BK) / 2,
    Frequency_Boer = sum(Frequency_Boer, na.rm = TRUE) / 2
  ) %>%
  ungroup()

#Round the frequencies to 3 decimal places
allele_freq <- allele_freq %>%
  mutate(
    Frequency_BK = round(Frequency_BK, 3),
    Frequency_Boer = round(Frequency_Boer, 3)
  )


#
# Calculate the total genotype frequencies for the entire population
geno_data_total <- geno_data %>%
  mutate(Total_Count = Count_BK + Count_Boer) %>%  # Sum counts across both breeds
  group_by(snpID) %>%
  mutate(Total_Individuals = sum(Total_Count, na.rm = TRUE)) %>%  # Calculate total individuals for each SNP
  ungroup() %>%
  mutate(Total_Frequency = Total_Count / Total_Individuals)  # Calculate frequency

# Summarize the results to get frequencies for each genotype
geno_data_total <- geno_data_total %>%
  select(snpID, alel, Total_Frequency) %>%
  arrange(snpID, alel)


