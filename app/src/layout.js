import React from 'react';
import { Outlet } from "react-router-dom";
import './layout.css'

function Layout() {
  return (
    <div>
      <nav className="navbar navbar-expand-lg navbar-light bg-light sticky-top">
        <a className="navbar-brand ontarget-logo" href="/">
          <img src="/ontarget_logo.svg" width="180" className="d-inline-block align-top" alt="" /></a>
        <button className="navbar-toggler" type="button" data-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
          <span className="navbar-toggler-icon"></span>
        </button>
        <div className="cnavbar" id="navbarNav">
          <ul className="navbar-nav">
            <li className="nav-item active">
              <a className="nav-link" href="/">Home</a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="/design">Design</a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="/pitx3example">PITX3 Examples</a>
            </li>
            <li className="nav-item">
              <a className="nav-link" href="/adora2example">ADORA2A Example</a>
            </li>
          </ul>
        </div>
      </nav>
      <Outlet />
    </div>
  );
}

export default Layout