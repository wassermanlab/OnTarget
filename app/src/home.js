import * as React from "react";
import { useNavigate } from "react-router-dom";

function Home() {

  let navigate = useNavigate(); 
  const routeChange = () =>{ 
    let path = `design`; 
    navigate(path);
  }
  

  return (
    <div className="container-fluid">
      <div className="row">
        <img src="/ontarget_infographic.svg" width="100%" className="d-inline-block align-top" alt="" />
      </div>
      <div className="row">
        <div className="centre">
          <img onClick={routeChange} src="/onTargetButton.svg" width="30%" className="d-inline-block align-top" alt="" />
        </div>
      </div>
    </div>
  );
}

export default Home