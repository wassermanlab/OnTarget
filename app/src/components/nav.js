import React from 'react';
import { Link } from "react-router-dom";
function Nav(props) {

  return (<React.Fragment>
    <nav className="navbar fixed-top navbar-expand-lg navbar-light bg-light">
      <Link className="navbar-brand" to="/">
      <img src={'/ontarget_logo.svg'} width="125"  alt="OnTarget" />
      </Link>
      <button className="navbar-toggler" type="button" data-toggle="collapse" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
        <span className="navbar-toggler-icon"></span>
      </button>
      <div className="collapse navbar-collapse" id="navbarNav">
        <ul className="navbar-nav">
          <li className="nav-item active">
            <Link className="nav-link" to="/about">
              About
            </Link>
          </li>
          <li className="nav-item active">
            <Link className="nav-link" to="/app">
              App
            </Link>
          </li>
        </ul>
      </div>

    </nav>



  </React.Fragment>)
}

export default Nav;