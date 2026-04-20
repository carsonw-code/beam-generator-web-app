# Beam Practice App

This is a Streamlit app I built to help statics students practice simply supported beam problems involving distributed loads, point loads, free moments, support reactions, shear diagrams, moment diagrams, and piecewise equations.

The goal of the app is to make beam-diagram practice faster and more interactive by generating new problems instantly and letting users reveal only the parts of the solution they want to check.

## Live app

https://beam-generator-web-app-mccunackxvmmeujrxljunr.streamlit.app/

## What it does

Each problem includes:
- one distributed load
- one point load
- one free moment

For every beam, the app can:
- generate a new random problem
- optionally reveal the support reactions
- reveal the shear and moment diagrams
- display the piecewise equations for each section

## Why I made it

I originally built this as a study tool for myself while learning statics, but I also wanted it to be useful for other students taking the course later on. The idea was to reduce the friction of practicing beam problems by making something quick to use, visual, and easy to check.

## How it works

For each new problem, the app:
1. generates a valid simply supported beam configuration
2. solves for the support reactions
3. splits the beam into sections based on loading changes
4. builds exact piecewise equations for \(V(x)\) and \(M(x)\)
5. plots the shear and moment diagrams directly from those equations

## Built with

- Python
- NumPy
- Matplotlib
- Streamlit

## Current scope

This version is designed for practice with:
- simply supported beams
- one distributed load
- one point load
- one free moment

  <img width="1333" height="741" alt="image" src="https://github.com/user-attachments/assets/ef0769cf-72bb-4861-9f36-3ead9e2d7cfc" />

<img width="1232" height="736" alt="image" src="https://github.com/user-attachments/assets/e6686aa6-afb2-4a2d-b962-602f47b6d85c" />


