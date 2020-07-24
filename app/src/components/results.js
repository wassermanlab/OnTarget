import React from 'react';
import NavButtons from './navButtons';
import IGV from './igv'

export default function Results(props) {

  return (
    <React.Fragment>
      <h1 className="text-center">Results</h1>
      <hr />
      <IGV genome={'hg19'} start={50000} end={500000} chrom={'chr1'} />
      <NavButtons next={props.next} back={props.back} page={props.page} />
    </React.Fragment>)
}