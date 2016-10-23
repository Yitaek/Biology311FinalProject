# Biology 311 Final Project: RosR Regulation
Contributors: Yitaek Hwang, Tosin Omofoye, Steven Pierre, Zhihao Zhu

## Abstract
Previous work by Sharma et al. and Tonner et al. have shown that VNG0258H or RosR is a transcription factor that regulates gene expression in Halobacterium Salinarum NRC-1 under extreme oxidative stress. The goal of this project was to verify the results of the aforementioned papers and delve deeper into identifying specific genes that RosR either activates or represses in response to oxidative stress. This allows to map the genes whose expression are directly affected by RosR. In comparing the gene expression ratio in the Δura3 parent and ΔVNG0258H mutant strains before and after exposing each subject to H2O2, we validated six expression profiled noted in Sharma et al. Then, we used k-means clustering to find a list of genes regulated by RosR using the same gene expression profile. This list, however, differed significantly from the list provided in Sharma et al., leading us to question our stringent criteria for ChIP-chip analysis. Next, we created a transcription factor network to find the relationship between the duration of H2O2 exposure to the genes that RosR regulates. The specific functions of these genes were then analyzed using the arCOG gene ontology analysis. Lastly, we compared our results with those of Tonner et al. to validate our analysis. In conclusion, the analysis showed that there are very few genes that have very high correlation to RosR binding at every stage of response to extreme peroxide levels, but many genes do have high correlation to RosR binding immediately after exposure to extreme levels of oxidation.

## File Structure
- Code: R and MATLAB code used to analyze the data
- Data: Gene expression and ChIP-chip data for Halobacterium Salinarum NRC-1
- Figures: MATLAB and R generated plots
- Research: Existing research on RosR
- Final Poster.pdf