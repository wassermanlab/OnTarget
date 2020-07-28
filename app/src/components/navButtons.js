import React from 'react';


export default function NavButtons(props) {

  return (
    <React.Fragment>
      {props.page === 3 ? <button type="button" onClick={props.submit} className="btn btn-nav">Submit</button> : <React.Fragment></React.Fragment>}
      {props.page !== 3 ? <button type="button" onClick={props.next} className="btn btn-nav">Next</button> : <React.Fragment></React.Fragment>}
      {props.page !== 1 ? <button type="button" onClick={props.back} className="btn btn-nav">Back</button> : <React.Fragment></React.Fragment>}
    </React.Fragment>)
}