// import React from 'react';
// import ReactDOM from 'react-dom/client';
import './index.css';
// import App from './App';
// import reportWebVitals from './reportWebVitals';

// const root = ReactDOM.createRoot(document.getElementById('root'));
// root.render(
//   <React.StrictMode>
//     <App />
//   </React.StrictMode>
// );

// // If you want to start measuring performance in your app, pass a function
// // to log results (for example: reportWebVitals(console.log))
// // or send to an analytics endpoint. Learn more: https://bit.ly/CRA-vitals
// reportWebVitals();

import { createRoot } from "react-dom/client";
// import Form from "@rjsf/core";
// import Form from '@rjsf/bootstrap-4';
// import App from './App';
import Form from '@rjsf/material-ui/v5';
// import { Button } from 'antd';

// import { withTheme } from '@rjsf/core';
// import { Theme as AntDTheme } from '@rjsf/antd';

// Make modifications to the theme with your own fields and widgets
//
// const Form = withTheme(AntDTheme);

// ReactDOM.render(
//   <React.StrictMode>
//     <App />
//   </React.StrictMode>,
//   document.getElementById('root')
// );

// const schema = {
//   title: "Todo",
//   type: "object",
//   required: ["title"],
//   properties: {
//     title: {type: "string", title: "Title", default: "A new task"},
//     done: {type: "boolean", title: "Done?", default: false}
//   }
// };

// const schema_txt = `{
//   "title": "Todo",
//   "type": "object",
//   "required": ["title"],
//   "properties": {
//     "title": {"type": "string", "title": "Title", "default": "A new task"},
//     "done": {"type": "boolean", "title": "Done?", "default": false}
//   }
// }`;

// const schema_txt = `{
//   "title": "A registration form",
//   "description": "A simple form example.",
//   "type": "object",
//   "required": [
//     "firstName",
//     "lastName"
//   ],
//   "properties": {
//     "firstName": {
//       "type": "string",
//       "title": "First name",
//       "default": "Chuck"
//     },
//     "lastName": {
//       "type": "string",
//       "title": "Last name"
//     },
//     "telephone": {
//       "type": "string",
//       "title": "Telephone",
//       "minLength": 10
//     },
//     "modes": {
//       "type": "array",
//       "items": { "$ref": "#/$defs/sub-array" },
//       "description": "Modes"
//     }
//   },
//   "$defs": {
//     "sub-array": {
//       "type": "object",
//       "properties": {
//         "field1": {
//           "type": "number",
//           "description": "field 1"
//         },
//         "field2": {
//           "type": "string",
//           "description": "field 2"
//         }
//       }
//     }
//   }
// }`;

const schema_txt = `{
  "title": "Schema dependencies",
  "description": "These samples are best viewed without live validation.",
  "type": "object",
  "properties": {
    "simple": {
      "src": "https://spacetelescope.github.io/understanding-json-schema/reference/object.html#dependencies",
      "title": "Simple",
      "type": "object",
      "properties": {
        "name": {
          "type": "string"
        },
        "credit_card": {
          "type": "number"
        }
      },
      "required": [
        "name"
      ],
      "dependencies": {
        "credit_card": {
          "properties": {
            "billing_address": {
              "type": "string"
            }
          },
          "required": [
            "billing_address"
          ]
        }
      }
    },
    "conditional": {
      "title": "Conditional",
      "$ref": "#/definitions/person"
    },
    "arrayOfConditionals": {
      "title": "Array of conditionals",
      "type": "array",
      "items": {
        "$ref": "#/definitions/person"
      }
    },
    "fixedArrayOfConditionals": {
      "title": "Fixed array of conditionals",
      "type": "array",
      "items": [
        {
          "title": "Primary person",
          "$ref": "#/definitions/person"
        }
      ],
      "additionalItems": {
        "title": "Additional person",
        "$ref": "#/definitions/person"
      }
    }
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
      "dependencies": {
        "Do you have any pets?": {
          "oneOf": [
            {
              "properties": {
                "Do you have any pets?": {
                  "enum": [
                    "No"
                  ]
                }
              }
            },
            {
              "properties": {
                "Do you have any pets?": {
                  "enum": [
                    "Yes: One"
                  ]
                },
                "How old is your pet?": {
                  "type": "number"
                }
              },
              "required": [
                "How old is your pet?"
              ]
            },
            {
              "properties": {
                "Do you have any pets?": {
                  "enum": [
                    "Yes: More than one"
                  ]
                },
                "Do you want to get rid of any?": {
                  "type": "boolean"
                }
              },
              "required": [
                "Do you want to get rid of any?"
              ]
            }
          ]
        }
      }
    }
  }
}`;

const schema = JSON.parse(schema_txt);

function onFormSubmit (event) {
  const json_str = JSON.stringify(event.formData);
  console.log("---Form submitted---");
  // console.log(event.formData);
  console.log(json_str);
};

function onFormChange (event) {
      console.log("---Form changed---");
      console.log(event.formData);
};

const log = (type) => console.log.bind(console, type);
const rootElement = document.getElementById("root");
const root = createRoot(rootElement);

root.render(
  <Form schema={schema}
        onChange={onFormChange}
        // onChange={log("changed")}
        // onSubmit={log("submitted")}
        onSubmit={onFormSubmit}
        onError={log("errors")} />
);

