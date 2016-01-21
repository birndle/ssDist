# ssDist

###Overview
 
 VEP plugin that writes distance to nearest donor and acceptor splice sites.
 
###Installation

 PERL5LIB=~/.vep/Plugins/:$PERL5LIB;
 mv ssDist.pm ~/.vep/Plugins

###Example

 perl variant_effect_predictor.pl -i variations.vcf --plugin ssDist
