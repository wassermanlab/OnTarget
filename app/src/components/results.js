import React from 'react';
import IGV from './igv'

export default function Results(props) {

  return (
    <React.Fragment>
      <IGV genome={'hg19'} start={50000} end={500000} chrom={'chr1'} />
      {/* TODO: add weights and recalculation of bigwig */}
    </React.Fragment>)
}