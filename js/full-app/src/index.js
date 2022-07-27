import './index.css';
import React from 'react';
import ReactDOM from 'react-dom';

// import Form from "@rjsf/core";
// import Form from '@rjsf/bootstrap-4';
import Form from '@rjsf/material-ui/v5';
// import Form from '@rjsf/antd';

const schema = {
  "title": "Schema dependencies",
  "description": "These samples are best viewed without live validation.",
  "type": "object",
  "properties": {
    "arrayOfConditionals": {
      "title": "Array of conditionals",
      "type": "array",
      "items": {
        "$ref": "#/definitions/person"
      }
    },
  },
  "definitions": {
    "person": {
      "title": "Person",
      "type": "object",
      "properties": {
        "Do you have any pets?": {
          "type": "string",
          "enum": [
            "No",
            "Yes: One",
            "Yes: More than one"
          ],
          "default": "No"
        }
      },
      "required": [
        "Do you have any pets?"
      ],
    }
  }
};

ReactDOM.render(<Form schema={schema} />,
        document.getElementById("root"));

