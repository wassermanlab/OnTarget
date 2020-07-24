import React, {useState} from 'react';
import Nav from '../components/nav'

export default function Home() {
    
    return (
        <React.Fragment>
            <Nav />
            <div className="navbar-fix">
                <img style={{ width: "100%" }} src={'/ontarget_infographic.svg'} alt="" />
            </div>
        </React.Fragment>)
}