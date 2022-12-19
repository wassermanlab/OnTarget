import React from 'react';

function NavButtons(props) {
  if (props.page === 1) {
    return (
    <div className='button-div'>
      <button type="button" className="btn btn-primary ontarget-button nav-button" onClick={props.next}>Next</button>
      </div>)
  }
  else if (props.page === 2) {
    return (
      <div className='button-div'>
        <button type="button" className="btn btn-primary ontarget-button nav-button" onClick={props.back}>Back</button> 
        <button type="button" className="btn btn-primary ontarget-button" onClick={props.getRegulatoryRegions}>Get Regulatory Regions</button>
        </div>)
  }
  return null;
}

export default NavButtons;