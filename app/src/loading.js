import React from 'react';

function Loading(props) {
  if (!props.loadedResources) {
  return (
    <React.Fragment>
      <p>Loading...</p>
    </React.Fragment>
  )
}
}
export default Loading;