# SailplanNg

This project was generated with [Angular CLI](https://github.com/angular/angular-cli) version 10.0.8.

## Development

Before you can build this project, you must install and configure the following dependencies on your machine:

1. [Node.js][]: We use Node to run a development web server and build the project.
   Depending on your system, you can install Node either from source or as a pre-packaged bundle.

After installing Node, you should be able to run the following command to install development tools.
You will only need to run this command when dependencies change in [package.json](package.json).

    npm ci

We use npm scripts and [Angular CLI][] as our build system.

## Development server

Run `npm start` for a dev server which will connect to the API assumed to be running locally at `http://localhost:8000`. Navigate to `http://localhost:3000/`. The app will automatically reload if you change any of the source files.

To connect to the Dev API instance run the following:

    npm run dev

To connect to the Production API instance run the following:

    npm run prod

To run all tooling prior to commiting source:

    npm run precommit

This runs linting, formatting, tests and does a production build. This is likely the command you'll want to run before you commit your changes.

## Code scaffolding

Run `ng generate component component-name` to generate a new component. You can also use `ng generate directive|pipe|service|class|guard|interface|enum|module`.

## Build

Run `ng build` to build the project. The build artifacts will be stored in the `dist/` directory. Use the `--prod` flag for a production build.

To run development build, run: `npm run build:dev`
To run production build, run: `npm run build:prod`

## Running unit tests

Run `npm run test` to execute the unit tests via [Jest](https://jestjs.io).

## Running end-to-end tests

Run `ng e2e` to execute the end-to-end tests via [Protractor](http://www.protractortest.org/).

## Code linting

`npm run lint` this runs linting on the project to report on programming errors/warnings
`npm run lint:fix` this runs linting and attempts to automatically fix a limited set of issues

## Code formatting

`npm run format` this runs prettier to indicate if the code is correctly formatted to our standard
`npm run format:fix` this runs prettier and automatically fixes the formatting on all source files

## Further help

To get more help on the Angular CLI use `ng help` or go check out the [Angular CLI README](https://github.com/angular/angular-cli/blob/master/README.md).
