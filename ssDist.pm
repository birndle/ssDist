=head1 CONTACT                                                                                                       

 Daniel Birnbaum <danpbirnbaum@gmail.com>
 
=cut

=head1 NAME

 ssDist

=head1 SYNOPSIS

 mv ssDist.pm ~/.vep/Plugins
 perl variant_effect_predictor.pl -i variations.vcf --plugin ssDist

=head1 DESCRIPTION

 A VEP plugin that writes the distance to nearest donor and acceptor splice sites for a given transcript.

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
    
    my @consequences = map { $_->SO_term } @{ $transcript_variation_allele->get_all_OverlapConsequences };
    my $genic_variant = !("upstream_gene_variant" ~~ @consequences || "downstream_gene_variant" ~~ @consequences);
    my $splice_lof = "splice_acceptor_variant" ~~ @consequences || "splice_donor_variant" ~~ @consequences;
    
    my ($donorDist, $acceptorDist);
    if ($genic_variant && !$splice_lof) {
        my ($dd, $ad) = get_dist_to_splice_sites($tv, $vf);        
        $donorDist = $dd if $dd ne 'NA';
        $acceptorDist = $ad if $ad ne 'NA';
    }
    return { donorDist => $donorDist, acceptorDist => $acceptorDist };
}


sub get_dist_to_splice_sites {
    my $tv = shift;
    my $vf = shift;

    my $slice = $vf->feature_Slice();
    my $tr = $tv->transcript;
    my $strand = $tr->strand();

    my ($dd, $da);
    if ($tv->exon_number) {
        my @exons = @{ $tr->get_all_Exons };
        my ($exon_num, $number_of_exons) = split /\//, ($tv->exon_number);
        my $exon = $exons[$exon_num - 1];
        $dd = ($strand == 1) ? $slice->{end} - $exon->{end} - 1 : $exon->{start} - $slice->{start} - 1;
        $da = ($strand == 1) ? $slice->{start} - $exon->{start} + 1 : $exon->{end} - $slice->{end} + 1;
        if ($exon_num == 1) {
            return (return_val($dd), 'NA');
        } elsif ($exon_num == $number_of_exons) {
            return ('NA', return_val($da));
        } else {
            return (return_val($dd), return_val($da));
        }
    } 
    elsif ($tv->intron_number) {
        my @introns = @{ $tr->get_all_Introns };
        my ($intron_num, $number_of_introns) = split /\//, ($tv->intron_number);
        my $intron = $introns[$intron_num - 1];
        $dd = ($strand == 1) ? $slice->{start} - $intron->{start} + 1 : $intron->{end} - $slice->{end} + 1;
        $da = ($strand == 1) ? $slice->{end} - $intron->{end} - 1 : $intron->{start} - $slice->{start} - 1;
        return (return_val($dd), return_val($da));
    }
    else {
        return check_all_junctions($tr, $slice, $strand);
    }
}

sub return_val {
    my $dist = shift;
    return ($dist < 0) ? $dist : '+' . $dist;
}


sub check_all_junctions {
    my ($tr, $slice, $strand, $intron) = @_[0..3];
    my $i = 0;
    my ($five_start, $five_end, $three_start, $three_end);
    foreach my $intron(@{$tr->get_all_Introns}) {
        if ($strand > 0) {
            ($five_start, $five_end) = ($intron->start - 3, $intron->start + 5);
            ($three_start, $three_end) = ($intron->end - 19, $intron->end + 3);
        } else {
            ($five_start, $five_end) = ($intron->end - 5, $intron->end + 3);
            ($three_start, $three_end) = ($intron->start - 3, $intron->start + 19);
        }
        if (overlap($slice->start, $slice->end, $five_start, $five_end)) {
            return ('INSERTION_AT_DONOR', 'NA');
        }
        if (overlap($slice->start, $slice->end, $three_start, $three_end)) {
            return ('NA', 'INSERTION_AT_ACCEPTOR');
        }
        $i++;
    }
    return ('NA', 'NA');
}


1;