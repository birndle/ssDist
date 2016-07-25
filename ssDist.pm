=head1 CONTACT                                                                                                       

 Daniel Birnbaum <danpbirnbaum@gmail.com>
 
=cut

=head1 NAME

 ssDist

=head1 SYNOPSIS

 mv ssDist.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin ssDist

=head1 DESCRIPTION

 A VEP plugin that computes the distance to the nearest donor and acceptor splice sites.

=cut

package ssDist;

use strict;
use warnings;

use Bio::EnsEMBL::Variation::Utils::BaseVepPlugin;
use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);
use Bio::Perl;
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);

sub get_header_info {
    return {
        acceptorDist => "Distance to nearest acceptor splice site",
        donorDist => "Distance to nearest donor splice site"
    };
}

sub feature_types {
    return ['Transcript'];
}

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);   
    return $self;
}

sub run {
    my ($self, $transcript_variation_allele) = @_;
    my $vf = $transcript_variation_allele->variation_feature;
    my $tv = $transcript_variation_allele->transcript_variation;
    my $allele = $transcript_variation_allele->allele_string();
    
    my @consequences = map { $_->SO_term } @{ $transcript_variation_allele->get_all_OverlapConsequences };
    my $genic_variant = !("upstream_gene_variant" ~~ @consequences || "downstream_gene_variant" ~~ @consequences);
    my $splice_lof = "splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences;
    my $indel = $allele =~ "-";
    if ($genic_variant && !($splice_lof && $indel)) {
        my ($dd, $ad) = get_dist_to_splice_sites($tv, $vf);        
        return { donorDist => $dd, acceptorDist => $ad };
    }
    return {};
}


sub get_dist_to_splice_sites {
    my ($tv, $vf) = @_[0..1];

    my $slice = $vf->feature_Slice();
    my $tr = $tv->transcript;
    my $strand = $tr->strand();

    my ($dd, $da);
    if ($tv->exon_number) {
        my @exons = @{ $tr->get_all_Exons };
        my ($exon_num, $number_of_exons) = split /\//, ($tv->exon_number);
        my $exon = $exons[$exon_num - 1];
        
        if ($strand == 1) {
            $dd = $slice->{end} - $exon->{end} - 1;
            $da = $slice->{start} - $exon->{start} + 1;
        } else {
            $dd = $exon->{start} - $slice->{start} - 1;
            $da = $exon->{end} - $slice->{end} + 1;
        }

        if ($exon_num == 1) {
            $da = 'FIRST_EXON';
        } elsif ($exon_num == $number_of_exons) {
            $dd = 'LAST_EXON';
        }

    } elsif ($tv->intron_number) {
        my @introns = @{ $tr->get_all_Introns };
        my ($intron_num, $number_of_introns) = split /\//, ($tv->intron_number);
        my $intron = $introns[$intron_num - 1];

        if ($strand == 1) {
            $dd = $slice->{start} - $intron->{start} + 1;
            $da = $slice->{end} - $intron->{end} - 1;
        } else {
            $dd = $intron->{end} - $slice->{end} + 1;
            $da = $intron->{start} - $slice->{start} - 1;
        }
    }
    
    else { # edge case: insertion occurring right at the splice junction
        ($dd, $da) = check_for_insertion_at_junction($tr, $slice);
    }
    return ($dd, $da);
}


sub check_for_insertion_at_junction {
    my ($tr, $slice) = @_[0..1];
    my $strand = $tr->strand();
    my ($five_start, $five_end, $three_start, $three_end);
    my @exons = @{ $tr->get_all_Exons };
    my @introns = @{ $tr->get_all_Introns };
    my $num_introns = scalar @introns;
    for (my $i=0; $i < $num_introns; $i++)  {
        my $intron = $introns[$i];
        if ($strand > 0) {
            ($five_start, $five_end) = ($intron->start - 3, $intron->start + 5);
            ($three_start, $three_end) = ($intron->end - 19, $intron->end + 3);
        } else {
            ($five_start, $five_end) = ($intron->end - 5, $intron->end + 3);
            ($three_start, $three_end) = ($intron->start - 3, $intron->start + 19);
        }
        if (overlap($slice->start, $slice->end, $five_start, $five_end)) {
            my $exon = $exons[$i];
            my $exon_length = $exon->end - $exon->start + 1;
            my $da = ($i == 0) ? 'FIRST_EXON' : $exon_length;
            return (0, $da);
        }
        if (overlap($slice->start, $slice->end, $three_start, $three_end)) {
            my $exon = $exons[$i + 1];
            my $exon_length = $exon->end - $exon->start + 1;
            my $dd = ($i == $num_introns - 1) ? 'LAST_EXON' : -$exon_length;
            return ($dd, 0);
        }
    }
    return ('') x 2;
}


1;