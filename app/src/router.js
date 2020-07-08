import React from 'react';
import {BrowserRouter as Router, Route} from "react-router-dom";
import Home from "./containers/Home";
import About from "./containers/About";
import App from "./containers/App";

function AppRouter() {
  return (
    <Router>
        <Route path="/" exact component={Home} />
        <Route path="/about/" component={About} />
        <Route path="/app" component={App} />
    </Router>
  );
}

export default AppRouter;