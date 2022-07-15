import { createRoot } from "react-dom/client";
// import Form from "@rjsf/core";
import Form from '@rjsf/bootstrap-4';
// import App from './App';

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

const schema_txt = `{
  "title": "A registration form",
  "description": "A simple form example.",
  "type": "object",
  "required": [
    "firstName",
    "lastName"
  ],
  "properties": {
    "firstName": {
      "type": "string",
      "title": "First name",
      "default": "Chuck"
    },
    "lastName": {
      "type": "string",
      "title": "Last name"
    },
    "telephone": {
      "type": "string",
      "title": "Telephone",
      "minLength": 10
    },
    "modes": {
      "type": "array",
      "items": { "$ref": "#/$defs/sub-array" },
      "description": "Modes"
    }
  },
  "$defs": {
    "sub-array": {
      "type": "object",
      "properties": {
        "field1": {
          "type": "number",
          "description": "field 1"
        },
        "field2": {
          "type": "string",
          "description": "field 2"
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
}

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

