import { createRoot } from "react-dom/client";
import Form from "@rjsf/core";
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

const schema_txt = `{
  "title": "Todo",
  "type": "object",
  "required": ["title"],
  "properties": {
    "title": {"type": "string", "title": "Title", "default": "A new task"},
    "done": {"type": "boolean", "title": "Done?", "default": false}
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
        // onChange={onFormChange}
        // onChange={log("changed")}
        // onSubmit={log("submitted")}
        onSubmit={onFormSubmit}
        onError={log("errors")} />
);

