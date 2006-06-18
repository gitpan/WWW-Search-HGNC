package WWW::Search::HGNC;
use base 'WWW::Search';

use warnings;
use strict;

use WWW::SearchResult;

use LWP::UserAgent;
use Text::RecordParser::Tab;
use Text::ParseWords;

our $VERSION = '0.01';

=head1 NAME

WWW::Search::HGNC - Access HGNC's database of proteins

=head1 SYNOPSIS

  use WWW::Search;
  my $search = new WWW::Search('HGNC');

  my @hgnc_ids = [ 9092, 12043 ];
  $search->native_query( \@hgnc_ids );

  while( my $prot = $search->next_result ) {
    printf "Approved symbol: %s\n", $prot->{approved_symbol};
    printf "Approved name: %s\n", $prot->{approved_name};
    printf "HGNC ID: %s\n", $prot->{hgnc_id};
    # ... etc.
  }

=head1 DESCRIPTION

The HUGO Gene Nomenclature Committee (HGNC) maintains a database of
human proteins (L<http://www.gene.ucl.ac.uk/nomenclature/>). This
module provides access to protein information via the L<WWW::Search>
interface.

=head1 RESULT FIELDS

The results returned by this module are L<WWW::SearchResult> objects
containing the following fields.

=head2 accession_numbers

  @values = @{ $prot->accession_numbers };

Corresponds to the 'Accession Numbers' HGNC field.

=head2 aliases

  @values = @{ $prot->aliases };

Corresponds to the 'Aliases' HGNC field.

=head2 approved_name

  $value = $prot->approved_name;

Corresponds to the 'Approved Name' HGNC field.

=head2 approved_symbol

  $value = $prot->approved_symbol;

Corresponds to the 'Approved Symbol' HGNC field.

=head2 chromosome

  $value = $prot->chromosome;

Corresponds to the 'Chromosome' HGNC field.

=head2 date_approved

  $value = $prot->date_approved;

Corresponds to the 'Date Approved' HGNC field.

=head2 date_modified

  $value = $prot->date_modified;

Corresponds to the 'Date Modified' HGNC field.

=head2 date_name_changed

  $value = $prot->date_name_changed;

Corresponds to the 'Date Name Changed' HGNC field.

=head2 entrez_gene_ids

  @values = @{ $prot->entrez_gene_ids };

Corresponds to the 'Entrez Gene ID' HGNC field.

=head2 mapped_entrez_gene_id

  $value = $prot->mapped_entrez_gene_id;

Corresponds to the 'Entrez Gene ID (mapped data)' HGNC field.

=head2 enzyme_ids

  @values = @{ $prot->enzyme_ids };

Corresponds to the 'Enzyme IDs' HGNC field.

=head2 mapped_gdb_id, gdb_id

  $value = $prot->mapped_gdb_id;
  $value = $prot->gdb_id;

Corresponds to the 'GDB ID (mapped data)' HGNC field.

=head2 gene_family_names

  @values = @{ $prot->gene_family_names };

Corresponds to the 'Gene Family Name' HGNC field.

=head2 hgnc_id

  $value = $prot->hgnc_id;

Corresponds to the 'HGNC ID' HGNC field.

=head2 locus_type

  $value = $prot->locus_type;

Corresponds to the 'Locus Type' HGNC field.

=head2 mgd_id

  $value = $prot->mgd_id;

Corresponds to the 'MGD ID' HGNC field.

=head2 misc_ids

  @values = @{ $prot->misc_ids };

Corresponds to the 'Misc IDs' HGNC field.

=head2 mapped_omim_id, omim_id

  $value = $prot->mapped_omim_id;
  $value = $prot->omim_id;

Corresponds to the 'OMIM ID (mapped data)' HGNC field.

=head2 previous_names

  $value = $prot->previous_names;

Corresponds to the 'Previous Names' HGNC field.

=head2 previous_symbols

  @values = @{ $prot->previous_symbols };

Corresponds to the 'Previous Symbols' HGNC field.

=head2 pubmed_ids

  @values = @{ $prot->pubmed_ids };

Corresponds to the 'Pubmed IDs' HGNC field.

=head2 mapped_refseq_id

  $value = $prot->mapped_refseq_id;

Corresponds to the 'RefSeq (mapped data)' HGNC field.

=head2 refseq_ids

  @values = @{ $prot->refseq_ids };

Corresponds to the 'RefSeq IDs' HGNC field.

=head2 status

  $value = $prot->status;

Corresponds to the 'Status' HGNC field.

=head2 mapped_uniprot_id, uniprot_id

  $value = $prot->mapped_uniprot_id;
  $value = $prot->uniprot_id;

Corresponds to the 'UniProt ID (mapped data)' HGNC field.

=cut

sub native_setup_search {
  my( $self, $query ) = @_;
  my @ids = ref $query eq 'ARRAY' ? @$query : ( $query );
  $self->{_ids} = \@ids;
  $self->{_idx} = 0;
  $self->user_agent(1);
}

=head2 native_retrieve_some

Fetches protein data from the Hugo Nomenclature Committee's database.

=cut

sub native_retrieve_some {
  my $self = shift;
  my $id = $self->{_ids}->[$self->{_idx}++] or return;
  my $url = $self->_url($id);

  my $data = $self->_fetch_data($url) or return;
  my $tp = new Text::RecordParser::Tab( data => $data, trim => 1 ) ;
  $tp->bind_header();

  my $rec = $tp->fetchrow_hashref;

  my %data = ( );
  my $fields = $self->_fields;
  while( my( $label, $field ) = each %$fields ) {
    my $value = $rec->{$label};
    $value = [ quotewords( '\s*,\s*', 0, $value ) ] if $field->{comma_delim};

    $data{$_} = $value foreach @{ $field->{fields} };
  }

  return undef unless $data{hgnc_id};

  $self->approximate_result_count(1);

  my $hit = new WWW::SearchResult();
  $hit->{$_} = $data{$_} for keys %data;
  $hit->url( $url );
  push @{$self->{cache}}, $hit;

  return 1;
}

sub _fetch_data {
  my( $self, $url ) = @_;
  my $res = $self->user_agent->get($url);
  return unless $res->is_success;
  return $res->content;
}

sub _fields { {
  'OMIM ID (mapped data)' => { fields => [ 'mapped_omim_id', 'omim_id' ] },
  'Date Name Changed' => { fields => [ 'date_name_changed' ] },
  'Previous Names' => { fields => [ 'previous_names' ] },
  'HGNC ID' => { fields => [ 'hgnc_id' ] },
  'Date Approved' => { fields => [ 'date_approved' ] },
  'Chromosome' => { fields => [ 'chromosome' ] },
  'Status' => { fields => [ 'status' ] },
  'Misc IDs' => { fields => [ 'misc_ids' ], comma_delim => 1 },
  'Approved Symbol' => { fields => [ 'approved_symbol' ] },
  'UniProt ID (mapped data)' => { fields => [ 'mapped_uniprot_id', 'uniprot_id' ] },
  'Accession Numbers' => { fields => [ 'accession_numbers' ], comma_delim => 1 },
  'Aliases' => { fields => [ 'aliases' ], comma_delim => 1 },
  'GDB ID (mapped data)' => { fields => [ 'mapped_gdb_id', 'gdb_id' ] },
  'Enzyme IDs' => { fields => [ 'enzyme_ids' ], comma_delim => 1 },
  'RefSeq IDs' => { fields => [ 'refseq_ids' ], comma_delim => 1 },
  'Approved Name' => { fields => [ 'approved_name' ] },
  'Previous Symbols' => { fields => [ 'previous_symbols' ], comma_delim => 1 },
  'Entrez Gene ID (mapped data)' => { fields => [ 'mapped_entrez_gene_id' ] },
  'MGD ID' => { fields => [ 'mgd_id' ] },
  'Gene Family Name' => { fields => [ 'gene_family_names' ], comma_delim => 1 },
  'Date Modified' => { fields => [ 'date_modified' ] },
  'Pubmed IDs' => { fields => [ 'pubmed_ids' ], comma_delim => 1 },
  'Locus Type' => { fields => [ 'locus_type' ] },
  'RefSeq (mapped data)' => { fields => [ 'mapped_refseq_id' ] },
  'Entrez Gene ID' => { fields => [ 'entrez_gene_ids' ], comma_delim => 1 },
} }

sub _url {
  my( $self, $id ) = @_;
  return 'http://www.gene.ucl.ac.uk/cgi-bin/nomenclature/gdlw.pl?title=Genew+output+data&hgnc_dbtag=on&col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_locus_type&col=gd_prev_sym&col=gd_prev_name&col=gd_aliases&col=gd_pub_chrom_map&col=gd_date2app_or_res&col=gd_date_mod&col=gd_date_name_change&col=gd_pub_acc_ids&col=gd_enz_ids&col=gd_pub_eg_id&col=gd_mgd_id&col=gd_other_ids&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=gd_gene_fam_name&col=md_gdb_id&col=md_eg_id&col=md_mim_id&col=md_refseq_id&col=md_prot_id&status=Approved&status=Approved+Non-Human&status=Entry+Withdrawn&status_opt=3&=on&where=gd_hgnc_id%3D'.$id.'&order_by=gd_app_sym_sort&limit=1&format=text&submit=submit&.cgifields=&.cgifields=status&.cgifields=chr&.cgifields=hgnc_dbtag';
}

=head1 AUTHOR

David Iberri, C<< <diberri at cpan.org> >>

=head1 BUGS

Please report any bugs or feature requests to
C<bug-www-hgnc at rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=WWW-Search-HGNC>.
I will be notified, and then you'll automatically be notified of progress on
your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc WWW::Search::HGNC

You can also look for information at:

=over 4

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/WWW-Search-HGNC>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/WWW-Search-HGNC>

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=WWW-Search-HGNC>

=item * Search CPAN

L<http://search.cpan.org/dist/WWW-Search-HGNC>

=back

=head1 COPYRIGHT & LICENSE

Copyright 2006 David Iberri, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.

=cut

1;
