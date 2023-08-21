document.getElementById('main-content__layout--console__form').onsubmit =
  function (event) {
    event.preventDefault();
    let cmd = document
      .getElementById('main-content__layout--console__form--cmd')
      .value.trim();
    console_process_cmd(cmd);
    document.getElementById('main-content__layout--console__form--cmd').value =
      '';
  };

// SCROLL TO BOTTOM
function console_scroll_to_bottom() {
  let outputs = document.querySelector(
    '.main-content__layout--console__outputs'
  );
  outputs.scrollTop = outputs.scrollHeight;
}

// ADD OUTPPUTS
function console_add_output(text) {
  document.querySelector(
    '.main-content__layout--console__outputs'
  ).innerHTML += `<pre class="console_line">${text}</pre><br><br>`;
  console_scroll_to_bottom();
}

function console_process_cmd(cmd) {
  let text;

  // HELP
  if (cmd == 'help') {
    text = [
      'CMD      DESCRIPTION',
      '----------------------------------',
      'help      Display this message.',
      'clear     Clear the console board.',
      'setup     Display setup configuration.',
      'solution  Show solution settings.',
      'results   Display numerical solutions.',
      'all       Show setup, solution and results.',
      'reset     Clear memory cache.',
      'about     Brief information.',
    ].join('<br>');
    console_add_output(text);
    return;
  }

  // CLEAR
  else if (cmd == 'clear') {
    document.querySelector(
      '.main-content__layout--console__outputs'
    ).innerHTML = '';
    return;
  }

  // SETUP
  else if (cmd == 'setup') {
    // MATERIAL PARAMETERS
    text = 'ZONES:<br>';
    text += '------<br>';
    if (ZN == 0) text += 'There are no zones specified yet.<br>';
    else {
      text += 'id    st       ss<br>';
      for (let i = 0; i < ZON.length; i++) {
        if (i < 10)
          text += `${i + 1}     ${ZON[i].st.toPrecision(4)}    ${ZON[
            i
          ].ss.toPrecision(4)}<br>`;
        else
          text += `${i + 1}    ${ZON[i].st.toPrecision(4)}    ${ZON[
            i
          ].ss.toPrecision(4)}<br>`;
      }
    }

    // REGIONS
    text += '<br>';
    text += 'REGIONS IN X:<br>';
    text += '-------------<br>';
    if (XREG == 0) text += 'There are no regions in x specified yet.<br>';
    else {
      text += 'id    len      nc<br>';
      for (let i = 0; i < XDOM.length; i++) {
        if (i < 10)
          text += `${i + 1}     ${XDOM[i].len.toPrecision(4)}    ${
            XDOM[i].nc
          }<br>`;
        else
          text += `${i + 1}    ${XDOM[i].len.toPrecision(4)}    ${
            XDOM[i].nc
          }<br>`;
      }
    }
    text += '<br>';
    text += 'REGIONS IN Y:<br>';
    text += '-------------<br>';
    if (YREG == 0) text += 'There are no regions in y specified yet.<br>';
    else {
      text += 'id    len      nc<br>';
      for (let i = 0; i < YDOM.length; i++) {
        if (i < 10)
          text += `${i + 1}     ${YDOM[i].len.toPrecision(4)}    ${
            YDOM[i].nc
          }<br>`;
        else
          text += `${i + 1}    ${YDOM[i].len.toPrecision(4)}    ${
            YDOM[i].nc
          }<br>`;
      }
    }

    // ZONE MAPPING
    text += '<br>';
    text += 'ZONE MAPPING (ZMAP):<br>';
    text += '--------------------<br>';
    if (ZN == 0 || XREG == 0 || YREG == 0) {
      text += "ZMAP can't be computed due to:<br>";
      text += '- There are no zones specified yet.<br>';
      text += '- There are no regions in x specified yet.<br>';
      text += '- There are no regions in y specified yet.<br>';
    } else {
      if (ZMAP == undefined) setup_zone_map_create();
      for (let j = YREG; j >= 0; j--) {
        for (let i = 0; i <= XREG; i++) {
          if (j == 0 && i == 0) text += 'RY\\RX';
          else if (j == 0) text += `  ${i}`;
          else if (i == 0) text += `${j}`;
          else if (i == 1 && j < 10) text += `      ${ZMAP[j - 1][i - 1] + 1}`;
          else if (i == 1 && j > 10) text += `     ${ZMAP[j - 1][i - 1] + 1}`;
          else text += `  ${ZMAP[j - 1][i - 1] + 1}`;
        }
        text += '<br>';
      }
    }

    // SOURCE MAPPING
    text += '<br>';
    text += 'SOURCE MAPPING (QMAP):<br>';
    text += '----------------------<br>';
    if (XREG == 0 || YREG == 0) {
      text += "QMAP can't be computed due to:<br>";
      text += '- There are no regions in x specified yet.<br>';
      text += '- There are no regions in y specified yet.<br>';
    } else {
      if (QMAP == undefined) setup_source_map_create();
      for (let j = YREG; j >= 0; j--) {
        for (let i = 0; i <= XREG; i++) {
          if (j == 0 && i == 0) text += 'RY\\RX';
          else if (j == 0) text += `      ${i}`;
          else if (i == 0) text += `${j}`;
          else if (i == 1 && j < 10)
            text += `      ${QMAP[j - 1][i - 1].toPrecision(4)}`;
          else if (i == 1 && j > 10)
            text += `     ${QMAP[j - 1][i - 1].toPrecision(4)}`;
          else text += `  ${QMAP[j - 1][i - 1].toPrecision(4)}`;
        }
        text += '<br>';
      }
    }

    // BOUNDARY CONDITIONS
    text += '<br>';
    text += 'BOUNDARY CONDITIONS:<br>';
    text += '--------------------<br>';

    if (BC['left'] == 0) text += 'Left:   Vacuum<br>';
    else if (BC['left'] == -1) text += 'Left:   Reflective<br>';
    else text += `Left:   Prescribed (value = ${BC['left']})<br>`;

    if (BC['right'] == 0) text += 'Right:  Vacuum<br>';
    else if (BC['right'] == -1) text += 'Right:  Reflective<br>';
    else text += `Right:  Prescribed (value = ${BC['right']})<br>`;

    if (BC['top'] == 0) text += 'Top:    Vacuum<br>';
    else if (BC['top'] == -1) text += 'Top:    Reflective<br>';
    else text += `Top:    Prescribed (value = ${BC['top']})<br>`;

    if (BC['bottom'] == 0) text += 'Bottom: Vacuum<br>';
    else if (BC['bottom'] == -1) text += 'Bottom: Reflective<br>';
    else text += `Bottom: Prescribed (value = ${BC['bottom']})<br>`;

    console_add_output(text);
    return;
  }

  // SOLUTION
  else if (cmd == 'solution') {
    // SELECT METHOD
    let options = document.querySelectorAll("input[name='scheme']");
    options.forEach((option) => {
      if (option.checked) {
        selected = option.id;
        METH = selected;
      }
    });

    // QUADRATURE
    text = 'QUADRATURE:<br>';
    text += '----------<br>';
    if (QUAD == undefined) text += 'To define<br><br>';
    else {
      let M = QUAD['miu'].length;
      text += ' MIU       THETA    W<br>';
      for (let m = 0; m < M; m++) {
        if (QUAD['miu'][m] > 0) {
          if (QUAD['theta'][m] > 0) {
            text += ` ${QUAD['miu'][m].toPrecision(5)}   ${QUAD['theta'][
              m
            ].toPrecision(5)}  ${QUAD['w'][m].toPrecision(5)}<br>`;
          } else {
            text += ` ${QUAD['miu'][m].toPrecision(5)}  ${QUAD['theta'][
              m
            ].toPrecision(5)}  ${QUAD['w'][m].toPrecision(5)}<br>`;
          }
        } else {
          if (QUAD['theta'][m] > 0) {
            text += `${QUAD['miu'][m].toPrecision(5)}   ${QUAD['theta'][
              m
            ].toPrecision(5)}  ${QUAD['w'][m].toPrecision(5)}<br>`;
          } else {
            text += `${QUAD['miu'][m].toPrecision(5)}  ${QUAD['theta'][
              m
            ].toPrecision(5)}  ${QUAD['w'][m].toPrecision(5)}<br>`;
          }
        }
      }
      text += '<br>';
    }

    // TOLERANCE
    text += `TOLERANCE = ${TOL}<br><br>`;

    // METHOD
    if (METH == 'dd_method')
      text += 'NUMERICAL SCHEME: <br>Diamond Difference (DD) Method<br>';
    else if (METH == 'step_method')
      text += 'NUMERICAL SCHEME: <br>Step (STP) Method<br>';
    else if (METH == 'ld_method')
      text += 'NUMERICAL SCHEME: <br>Linear Discontinuous (LD) Method<br>';
    else if (METH == 'rm_cn_method')
      text +=
        'NUMERICAL SCHEME: <br>Response Matrix - Constant Nodal (RM-CN) Method<br>';
    else if (METH == 'rm_lln_method')
      text +=
        'NUMERICAL SCHEME: <br>Response Matrix - Linear Linear Nodal (RM-LLN) Method<br>';

    console_add_output(text);
  }

  // RESULTS
  else if (cmd == 'results') {
    // VALIDATION
    if (results_validator() == false) {
      text = 'No results detected<br><br>';
    } else {
      // ITERATIONS
      text = `ITER: ${ITER}<br>`;
      text += '------------<br><br>';

      // CPU
      text += `CPU: ${CPU.toExponential(4)}<br>`;
      text += '------------<br><br>';

      // SCALAR FLUX
      text += 'SCALAR FLUX:<br>';
      text += '------------<br>';
      let flux_r = results_report_data('scalar_flux');
      for (let j = YREG; j >= 0; j--) {
        for (let i = 0; i <= XREG; i++) {
          if (j == 0 && i == 0) text += 'RY\\RX';
          else if (j == 0 && i != 0) {
            if (i == 1) text += `  ${i}`;
            else text += `          ${i}`;
          } else if (i == 0) text += `${j}`;
          else if (i == 1 && j < 10) {
            if (flux_r[j - 1][i - 1] > 0.0) {
              text += `      ${flux_r[j - 1][i - 1].toExponential(4)}`;
            } else {
              text += `     ${flux_r[j - 1][i - 1].toExponential(4)}`;
            }
          } else if (i == 1 && j > 10) {
            if (flux_r[j - 1][i - 1] > 0.0) {
              text += `      ${flux_r[j - 1][i - 1].toExponential(4)}`;
            } else {
              text += `     ${flux_r[j - 1][i - 1].toExponential(4)}`;
            }
          } else {
            if (flux_r[j - 1][i - 1] > 0.0)
              text += `  ${flux_r[j - 1][i - 1].toExponential(4)}`;
            else text += ` ${flux_r[j - 1][i - 1].toExponential(4)}`;
          }
        }
        text += '<br>';
      }

      // LEAKAGE
      text += '<br>';
      text += 'LEAKAGE:<br>';
      text += '--------<br>';
      text += 'LEFT       BOTTOM     RIGHT      TOP<br>';
      let leak = results_report_data('leakage');
      // LEFT
      if (BC['left'] == -1) text += `-          `;
      else if (leak[0] > 0.0) text += `${leak[0].toExponential(4)}  `;
      else text += `${leak[0].toExponential(4)} `;
      // BOTTOM
      if (BC['bottom'] == -1) text += `-          `;
      else if (leak[1] > 0.0) text += `${leak[1].toExponential(4)}  `;
      else text += `${leak[1].toExponential(4)} `;
      // RIGHT
      if (BC['right'] == -1) text += `-          `;
      else if (leak[2] > 0.0) text += `${leak[2].toExponential(4)}  `;
      else text += `${leak[2].toExponential(4)} `;
      // BOTTOM
      if (BC['top'] == -1) text += `-          `;
      else if (leak[3] > 0.0) text += `${leak[3].toExponential(4)}  `;
      else text += `${leak[3].toExponential(4)} `;
    }
    console_add_output(text);
  }

  // ALL
  else if (cmd == 'all') {
    console_process_cmd('clear');
    console_process_cmd('setup');
    console_process_cmd('solution');
    console_process_cmd('results');
  } else if (cmd == 'reset') {
    // SETUP
    ZON = 0;
    ZON = [];
    XREG = 0;
    XDOM = [];
    YREG = 0;
    YDOM = [];
    ZMAP = undefined;
    QMAP = undefined;
    BC = {left: 0, right: 0, top: 0, bottom: 0};

    // SOLUTION
    QUAD = undefined;
    TOL = 0.00001;
    METH = 'dd_method';

    // RESULTS
    results_reset();

    if (document.getElementById('canvas_geom'))
      document.getElementById('canvas_geom').remove();
    if (document.getElementById('canvas_quad'))
      document.getElementById('canvas_quad').remove();
    if (document.getElementById('results__picture'))
      document.getElementById('results__picture').remove();
    document.getElementById('initial_img').style.display = 'block';
  }

  // ABOUT
  else if (cmd == 'about') {
    text = [
      'NEUTRON TRANSPORT SIMULATOR (NTS) v2.0',
      '--------------------------------------',
      "VIDEO DEMO: <a href='https://youtu.be/AtY52V7GDaI' target='_blank'>https://youtu.be/AtY52V7GDaI</a>",
    ].join('<br>');
    console_add_output(text);
    return;
  } else if (cmd == 'calculation_done') {
    text = 'Program terminated with status code 00:<br>';
    text += '---------------------------------------<br>';
    text += 'Calculation done';

    console_add_output(text);
    return;
  }

  // EMPTY
  else if (cmd == '') {
    return;
  }

  // OTHER
  else {
    text = 'Command not found.<br>';
    console_add_output(text);
    return;
  }
}

function console_process_error(error) {
  let text;

  // GEOMETRY ERROR
  if (error == 'geometry') {
    text = "Geometry can't be displayed due to:<br>";
    text += '- There are no zones specified yet.<br>';
    text += '- There are no regions in x specified yet.<br>';
    text += '- There are no regions in y specified yet.<br>';
    text += "- ZMAP isn't specified yet.<br>";
    text += "- QMAP isn't specified yet.<br>";

    console_add_output(text);
    return;
  }

  // CALCULATION ERROR DUE TO GEOMETRY
  else if (error == 'calculation') {
    text = "Calculation can't be done due to:<br>";
    text += '- There are no zones specified yet.<br>';
    text += '- There are no regions in x specified yet.<br>';
    text += '- There are no regions in y specified yet.<br>';
    text += "- ZMAP isn't specified yet.<br>";
    text += "- QMAP isn't specified yet.<br>";

    console_add_output(text);
    return;
  }

  // CALCULATION ERROR DUE TO QUADRATURE
  else if (error == 'quadrature') {
    text = "Calculation can't be done due to:<br>";
    text += '- Quadrature order not specified yet.<br>';

    console_add_output(text);
    return;
  }

  // CALCULATION ERROR DUE TO RESULTS
  else if (error == 'results') {
    text = 'Results problem due to:<br>';
    text += '- No results detected.<br>';
    text += '- Data incompatibility.<br>';

    console_add_output(text);
    return;
  }

  // ERROR INSIDE DE SERVER
  else if (error == '1') {
    text = 'Program terminated with status code 01:<br>';
    text += '---------------------------------------<br>';
    text += 'Program usage (server side)';

    console_add_output(text);
    return;
  } else if (error == '2') {
    text = 'Program terminated with status code 02:<br>';
    text += '---------------------------------------<br>';
    text += 'Incompatible problem configuration (client side)';

    console_add_output(text);
    return;
  } else if (error == '3') {
    text = 'Program terminated with status code 03:<br>';
    text += '---------------------------------------<br>';
    text += 'Memory allocation problem (server side)';

    console_add_output(text);
    return;
  } else if (error == '4') {
    text = 'Program terminated with status code 04:<br>';
    text += '---------------------------------------<br>';
    text += 'Singular response matrix (client side)';

    console_add_output(text);
    return;
  } else if (error == '5') {
    text = 'Program terminated with status code 05:<br>';
    text += '---------------------------------------<br>';
    text += 'Maximum iteration number reached (client side)';

    console_add_output(text);
    return;
  } else if (error == '6') {
    text = 'Program terminated with status code 06:<br>';
    text += '---------------------------------------<br>';
    text += 'Arithmetic precision problem (client side)';

    console_add_output(text);
    return;
  } else if (error == '7') {
    text = 'Program terminated with status code 07:<br>';
    text += '---------------------------------------<br>';
    text += 'Server problem detected';

    console_add_output(text);
    return;
  } else if (error == '8') {
    text = 'Program terminated with status code 08:<br>';
    text += '---------------------------------------<br>';
    text += 'Bad request (client side)';

    console_add_output(text);
    return;
  }

  return;
}
