/////////////////////////////////////////////////////////////////////// RESULTS SECTION => VARIABLES
/////////////////////////////////////////////////////////////////////////////////////////////////////
let MFLUX;
let MFLOW;
let XFLOW;
let YFLOW;
let ITER;
let CPU;

// RESET RESULTS
function results_reset() {
  MFLUX = undefined;
  MFLOW = undefined;
  XFLOW = undefined;
  YFLOW = undefined;
  ITER = undefined;
  CPU = undefined;
}

///////////////////////////////////////////////////////////////////// RESULTS SECTION => INITIALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////
function results_section_init() {
  // CLEAR UPLOAD RESULTS
  results_upload_file_reset();
}

////////////////////////////////////////////////////////////////////// RESULTS SECTION => UPLOAD FILE
/////////////////////////////////////////////////////////////////////////////////////////////////////

// FUNCTION FOR LOAD EXTERNAL FILE
function results_upload_file() {
  let input = document.getElementById('results-section__input--upload-file');
  let files = input.files;
  if (files.length == 0) return;
  const file = files[0];
  let reader = new FileReader();
  reader.onload = function (e) {
    let input = document.getElementById('results-section__input--upload-file');
    let err = document.getElementById(
      'results-section__input--upload-file__msg'
    );
    let data;
    try {
      data = JSON.parse(e.target.result);
    } catch (e) {
      input.classList.remove('is-valid');
      input.classList.remove('is-invalid');
      input.classList.add('is-invalid');
      err.className = 'invalid-feedback';
      err.innerHTML = 'Invalid file';
      return;
    }
    if (results_upload_file_validator(data)) {
      // SETUP
      ZN = data['ZN'];
      ZON = data['ZON'];
      XREG = data['XREG'];
      XDOM = data['XDOM'];
      YREG = data['YREG'];
      YDOM = data['YDOM'];
      ZMAP = data['ZMAP'];
      QMAP = data['QMAP'];
      BC = data['BC'];

      // SOLUTION
      QUAD = data['QUAD'];
      TOL = data['TOL'];
      METH = data['METH'];

      // RESULTS
      ITER = data['ITER'];
      CPU = data['CPU'];
      MFLUX = data['MFLUX'];
      MFLOW = data['MFLOW'];
      XFLOW = data['XFLOW'];
      YFLOW = data['YFLOW'];

      input.classList.remove('is-valid');
      input.classList.remove('is-invalid');
      input.classList.add('is-valid');
      err.className = 'valid-feedback';
      err.innerHTML = 'Upload Completed';
    } else {
      input.classList.remove('is-valid');
      input.classList.remove('is-invalid');
      input.classList.add('is-invalid');
      err.className = 'invalid-feedback';
      err.innerHTML = 'Invalid data';
    }
  };
  reader.readAsText(file);
}

// FUNCTION FOR RESET FILE INPUT
function results_upload_file_reset() {
  let input = document.getElementById('results-section__input--upload-file');
  let err = document.getElementById('results-section__input--upload-file__msg');
  input.value = '';
  input.classList.remove('is-valid');
  input.classList.remove('is-invalid');
  err.className = '';
  err.innerHTML = '';
}

// FUNCTION TO CHECK FOR INVALID DATAFILE
function results_upload_file_validator(data) {
  // ZONES
  if (data['ZN'] <= 0 || data['ZN'] == undefined) return false;
  if (data['ZON'] == undefined || data['ZON'].length != data['ZN'])
    return false;
  for (let z = 0; z < data['ZON'].length; z++) {
    if (data['ZON'][z].st <= 0 || data['ZON'][z].st == undefined) return false;
    if (
      data['ZON'][z].ss < 0 ||
      data['ZON'][z].ss == undefined ||
      data['ZON'][z].ss >= data['ZON'][z].st
    )
      return false;
  }

  // X REGIONS
  if (data['XREG'] <= 0 || data['XREG'] == undefined) return false;
  if (data['XDOM'] <= undefined || data['XDOM'].length != data['XREG'])
    return false;
  for (let i = 0; i < data['XREG'].length; i++) {
    if (data['XDOM'][i].len <= 0 || data['XDOM'][i].len == undefined)
      return false;
    if (
      data['XDOM'][i].nc < 0 ||
      data['XDOM'][i].nc == undefined ||
      !Number.isInteger(data['XDOM'][i].nc)
    )
      return false;
  }

  // Y REGIONS
  if (data['YREG'] <= 0 || data['YREG'] == undefined) return false;
  if (data['YDOM'] <= undefined || data['YDOM'].length != data['YREG'])
    return false;
  for (let i = 0; i < data['YREG'].length; i++) {
    if (data['YDOM'][i].len <= 0 || data['YDOM'][i].len == undefined)
      return false;
    if (
      data['YDOM'][i].nc < 0 ||
      data['YDOM'][i].nc == undefined ||
      !Number.isInteger(data['YDOM'][i].nc)
    )
      return false;
  }

  // ZMAP
  if (data['ZMAP'] == undefined || data['ZMAP'].length != data['YREG'])
    return false;
  if (data['ZMAP'][0] == undefined || data['ZMAP'][0].length != data['XREG'])
    return false;
  for (let i = 0; i < data['XREG'].length; i++) {
    for (let j = 0; j < data['YREG'].length; j++) {
      if (
        data['ZMAP'][j][i] < 0 ||
        data['ZMAP'][j][i] >= data['ZN'] ||
        !Number.isInteger(data['ZMAP'][j][i])
      )
        return false;
    }
  }

  // QMAP
  if (data['QMAP'] == undefined || data['QMAP'].length != data['YREG'])
    return false;
  if (data['QMAP'][0] == undefined || data['QMAP'][0].length != data['XREG'])
    return false;
  for (let i = 0; i < data['XREG'].length; i++) {
    for (let j = 0; j < data['YREG'].length; j++) {
      if (data['QMAP'][j][i] < 0) return false;
    }
  }

  // BC
  if (data['BC'] == undefined) return false;
  if (
    data['BC'].left == undefined ||
    (data['BC'].left < 0 && data['BC'].left != -1)
  )
    return false;
  if (
    data['BC'].right == undefined ||
    (data['BC'].right < 0 && data['BC'].right != -1)
  )
    return false;
  if (
    data['BC'].top == undefined ||
    (data['BC'].top < 0 && data['BC'].top != -1)
  )
    return false;
  if (
    data['BC'].bottom == undefined ||
    (data['BC'].bottom < 0 && data['BC'].bottom != -1)
  )
    return false;

  // QUADRATURE
  if (data['QUAD'] == undefined) return false;
  if (
    data['QUAD'].miu == undefined ||
    data['QUAD'].theta == undefined ||
    data['QUAD'].w == undefined
  )
    return false;
  let quad_len = data['QUAD'].miu.length;
  if (
    data['QUAD'].miu.length != quad_len ||
    data['QUAD'].theta.length != quad_len ||
    data['QUAD'].w.length != quad_len
  )
    return false;

  // TOLERANCE
  if (data['TOL'] == undefined || data['TOL'] >= 1 || data['TOL'] <= 0)
    return false;

  // METHOD
  if (data['METH'] == undefined) return false;

  // ITERATIONS
  if (data['ITER'] == undefined) return false;

  // CPU TIME
  if (data['CPU'] == undefined) return false;

  // SCALAR FLUX
  if (data['MFLUX'] == undefined) return false;
  let ntcx = 0,
    ntcy = 0;
  for (let rx = 0; rx < data['XREG']; rx++) ntcx = ntcx + data['XDOM'][rx].nc;
  for (let ry = 0; ry < data['YREG']; ry++) ntcy = ntcy + data['YDOM'][ry].nc;
  if (data['MFLUX'][0].length != ntcx || data['MFLUX'].length != ntcy)
    return false;

  // MEAN ANGULAR FLUX
  if (data['MFLOW'] == undefined) return false;
  if (
    data['MFLOW'].length != quad_len ||
    data['MFLOW'][0].length != ntcy ||
    data['MFLOW'][0][0].length != ntcx
  )
    return false;

  // ANGULAR FLUX IN X
  if (data['XFLOW'] == undefined) return false;
  if (
    data['XFLOW'].length != quad_len ||
    data['XFLOW'][0].length != ntcy ||
    data['XFLOW'][0][0].length != ntcx + 1
  )
    return false;

  // ANGULAR FLUX IN Y
  if (data['YFLOW'] == undefined) return false;
  if (
    data['YFLOW'].length != quad_len ||
    data['YFLOW'][0].length != ntcy + 1 ||
    data['YFLOW'][0][0].length != ntcx
  )
    return false;

  return true;
}

/////////////////////////////////////////////////////////////////////////// RESULTS SECTION => REPORT
/////////////////////////////////////////////////////////////////////////////////////////////////////
function results_validator() {
  // ZONES
  if (ZN <= 0 || ZN == undefined) return false;
  if (ZON == undefined || ZON.length != ZN) return false;
  for (let z = 0; z < ZON.length; z++) {
    if (ZON[z].st <= 0 || ZON[z].st == undefined) return false;
    if (ZON[z].ss < 0 || ZON[z].ss == undefined || ZON[z].ss >= ZON[z].st)
      return false;
  }

  // X REGIONS
  if (XREG <= 0 || XREG == undefined) return false;
  if (XDOM <= undefined || XDOM.length != XREG) return false;
  for (let i = 0; i < XREG; i++) {
    if (XDOM[i].len <= 0 || XDOM[i].len == undefined) return false;
    if (
      XDOM[i].nc < 0 ||
      XDOM[i].nc == undefined ||
      !Number.isInteger(XDOM[i].nc)
    )
      return false;
  }

  // Y REGIONS
  if (YREG <= 0 || YREG == undefined) return false;
  if (YDOM <= undefined || YDOM.length != YREG) return false;
  for (let i = 0; i < YREG; i++) {
    if (YDOM[i].len <= 0 || YDOM[i].len == undefined) return false;
    if (
      YDOM[i].nc < 0 ||
      YDOM[i].nc == undefined ||
      !Number.isInteger(YDOM[i].nc)
    )
      return false;
  }

  // ZMAP
  if (ZMAP == undefined || ZMAP.length != YREG) return false;
  if (ZMAP[0] == undefined || ZMAP[0].length != XREG) return false;
  for (let i = 0; i < XREG; i++) {
    for (let j = 0; j < YREG; j++) {
      if (ZMAP[j][i] < 0 || ZMAP[j][i] >= ZN || !Number.isInteger(ZMAP[j][i]))
        return false;
    }
  }

  // QMAP
  if (QMAP == undefined || QMAP.length != YREG) return false;
  if (QMAP[0] == undefined || QMAP[0].length != XREG) return false;
  for (let i = 0; i < XREG; i++) {
    for (let j = 0; j < YREG; j++) {
      if (QMAP[j][i] < 0) return false;
    }
  }

  // BC
  if (BC == undefined) return false;
  if (BC.left == undefined || (BC.left < 0 && BC.left != -1)) return false;
  if (BC.right == undefined || (BC.right < 0 && BC.right != -1)) return false;
  if (BC.top == undefined || (BC.top < 0 && BC.top != -1)) return false;
  if (BC.bottom == undefined || (BC.bottom < 0 && BC.bottom != -1))
    return false;

  // QUADRATURE
  if (QUAD == undefined) return false;
  if (QUAD.miu == undefined || QUAD.theta == undefined || QUAD.w == undefined)
    return false;
  let quad_len = QUAD.miu.length;
  if (
    QUAD.miu.length != quad_len ||
    QUAD.theta.length != quad_len ||
    QUAD.w.length != quad_len
  )
    return false;

  // TOLERANCE
  if (TOL == undefined || TOL >= 1 || TOL <= 0) return false;

  // METHOD
  if (METH == undefined) return false;

  // ITERATIONS
  if (ITER == undefined) return false;

  // CPU TIME
  if (CPU == undefined) return false;

  // SCALAR FLUX
  if (MFLUX == undefined) return false;
  let ntcx = 0,
    ntcy = 0;
  for (let rx = 0; rx < XREG; rx++) ntcx = ntcx + XDOM[rx].nc;
  for (let ry = 0; ry < YREG; ry++) ntcy = ntcy + YDOM[ry].nc;
  if (MFLUX[0].length != ntcx || MFLUX.length != ntcy) return false;

  // MEAN ANGULAR FLUX
  if (MFLOW == undefined) return false;
  if (
    MFLOW.length != quad_len ||
    MFLOW[0].length != ntcy ||
    MFLOW[0][0].length != ntcx
  )
    return false;

  // ANGULAR FLUX IN X
  if (XFLOW == undefined) return false;
  if (
    XFLOW.length != quad_len ||
    XFLOW[0].length != ntcy ||
    XFLOW[0][0].length != ntcx + 1
  )
    return false;

  // ANGULAR FLUX IN Y
  if (YFLOW == undefined) return false;
  if (
    YFLOW.length != quad_len ||
    YFLOW[0].length != ntcy + 1 ||
    YFLOW[0][0].length != ntcx
  )
    return false;

  return true;
}

// REPORT EVENT
document.getElementById('results-section__report--btn').onclick = function () {
  let options = document.querySelectorAll("input[name='report']");
  results_report_init(options);
};

// INITIALIZE REPORT
function results_report_init(options) {
  let msg = document.querySelector('#modal-report__msg');
  if (MFLUX == undefined) {
    msg.style.display = 'block';
    msg.innerHTML = 'There are no results loaded';
    return;
  } else if (results_validator() == false) {
    msg.style.display = 'block';
    msg.innerHTML = 'Data are not compatible';
    return;
  }

  let modal = document.getElementById('modal-report__content');
  modal.innerHTML = '';

  let iter = document.createElement('div');
  iter.className = 'my-2';
  iter.innerHTML = `<span class="font-weight-bold float-left">Number of Iterations: ${ITER}</span>`;
  modal.append(iter);

  let cpu = document.createElement('div');
  cpu.className = 'my-2';
  cpu.innerHTML = `<span class="font-weight-bold float-left">CPU Time (s): ${CPU.toExponential(
    5
  )}</span>`;
  modal.append(cpu);

  options.forEach((option) => {
    if (option.checked) {
      selected = option.id;

      // SCALAR FLUX PER REGION
      if (selected == 'results-section__report--scalar-flux') {
        let flux_r = results_report_data('scalar_flux');

        // TITLE
        let title = document.createElement('div');
        title.className = 'my-2';
        title.innerHTML =
          '<span class="font-weight-bold float-left">Scalar flux per region [n cm<sup>-2</sup> s<sup>-1</sup>]:</span>';
        modal.append(title);

        // TABLE
        let container = document.createElement('div');
        container.className = 'table-responsive mb-4';
        let table = document.createElement('table');
        table.className = 'table table-bordered table-sm';
        let tbody = document.createElement('tbody');
        for (let j = YREG; j >= 0; j--) {
          let tr = document.createElement('tr');
          for (let i = 0; i <= XREG; i++) {
            let td = document.createElement('td');
            if (i == 0 && j == 0) td.innerHTML = 'RY \\ RX';
            else if (j == 0) td.innerHTML = `${i}`;
            else if (i == 0) td.innerHTML = `${j}`;
            else td.innerHTML = `${flux_r[j - 1][i - 1].toExponential(5)}`;
            tr.append(td);
          }
          tbody.append(tr);
        }
        table.append(tbody);
        container.append(table);
        modal.append(container);
      }

      // ABSORPTION RATE PER REGION
      else if (selected == 'results-section__report--absorption-rate') {
        let abs_r = results_report_data('absorption_rate');

        // TITLE
        let title = document.createElement('div');
        title.className = 'my-2';
        title.innerHTML =
          '<span class="font-weight-bold float-left">Absorption rate per region [n s<sup>-1</sup>]:</span>';
        modal.append(title);

        // TABLE
        let container = document.createElement('div');
        container.className = 'table-responsive mb-4';
        let table = document.createElement('table');
        table.className = 'table table-bordered table-sm';
        let tbody = document.createElement('tbody');
        for (let j = YREG; j >= 0; j--) {
          let tr = document.createElement('tr');
          for (let i = 0; i <= XREG; i++) {
            let td = document.createElement('td');
            if (i == 0 && j == 0) td.innerHTML = 'RY \\ RX';
            else if (j == 0) td.innerHTML = `${i}`;
            else if (i == 0) td.innerHTML = `${j}`;
            else td.innerHTML = `${abs_r[j - 1][i - 1].toExponential(5)}`;
            tr.append(td);
          }
          tbody.append(tr);
        }
        table.append(tbody);
        container.append(table);
        modal.append(container);
      }

      // LEAKAGES
      else {
        let leak = results_report_data('leakage');

        // TITLE
        let title = document.createElement('div');
        title.className = 'my-2';
        title.innerHTML =
          '<span class="font-weight-bold float-left">Leakages at the boundaries [n cm<sup>-1</sup> s<sup>-1</sup>]:</span>';
        modal.append(title);

        // TABLE
        let container = document.createElement('div');
        container.className = 'table-responsive mb-4';
        let table = document.createElement('table');
        table.className = 'table table-bordered table-sm';
        let tbody = document.createElement('tbody');
        let trh = document.createElement('tr');
        trh.innerHTML =
          '<td>Left</td><td>Bottom</td><td>Right</td><td>Top</td>';
        tbody.append(trh);
        let tr = document.createElement('tr');
        let td0 = document.createElement('td');
        if (BC['left'] != -1) td0.innerHTML = leak[0].toExponential(5);
        else td0.innerHTML = '-';
        tr.append(td0);
        let td1 = document.createElement('td');
        if (BC['bottom'] != -1) td1.innerHTML = leak[1].toExponential(5);
        else td1.innerHTML = '-';
        tr.append(td1);
        let td2 = document.createElement('td');
        if (BC['right'] != -1) td2.innerHTML = leak[2].toExponential(5);
        else td2.innerHTML = '-';
        tr.append(td2);
        let td3 = document.createElement('td');
        if (BC['top'] != -1) td3.innerHTML = leak[3].toExponential(5);
        else td3.innerHTML = '-';
        tr.append(td3);
        tbody.append(tr);
        table.append(tbody);
        container.append(table);
        modal.append(container);
      }
    }
  });
}

// GENERATE REPORT DATA
function results_report_data(selected) {
  // SCALAR FLUX PER REGION
  if (selected == 'scalar_flux') {
    let i_b, len_x, nc_x, hi, j_b, len_y, nc_y, hj, area;
    let flux_r = [];
    for (let yr = 0; yr < YREG; yr++) {
      let aux = [];
      for (let xr = 0; xr < XREG; xr++) {
        aux.push(0);
      }
      flux_r.push(aux);
    }
    j_b = 0;
    for (let yr = 0; yr < YREG; yr++) {
      nc_y = YDOM[yr].nc;
      len_y = YDOM[yr].len;
      hj = len_y / nc_y;
      for (let j = 0; j < nc_y; j++) {
        i_b = 0;
        for (let xr = 0; xr < XREG; xr++) {
          nc_x = XDOM[xr].nc;
          len_x = XDOM[xr].len;
          hi = len_x / nc_x;
          area = len_x * len_y;
          for (let i = 0; i < nc_x; i++) {
            flux_r[yr][xr] =
              flux_r[yr][xr] + (hj * hi * MFLUX[j_b][i_b]) / area;
            i_b = i_b + 1;
          }
        }
        j_b = j_b + 1;
      }
    }
    return flux_r;
  }

  // ABSORPTION RATE PER REGION
  else if (selected == 'absorption_rate') {
    let i_b, len_x, nc_x, hi, j_b, len_y, nc_y, hj, z, sa;
    let abs_r = [];
    for (let yr = 0; yr < YREG; yr++) {
      let aux = [];
      for (let xr = 0; xr < XREG; xr++) {
        aux.push(0);
      }
      abs_r.push(aux);
    }
    j_b = 0;
    for (let yr = 0; yr < YREG; yr++) {
      nc_y = YDOM[yr].nc;
      len_y = YDOM[yr].len;
      hj = len_y / nc_y;
      for (let j = 0; j < nc_y; j++) {
        i_b = 0;
        for (let xr = 0; xr < XREG; xr++) {
          nc_x = XDOM[xr].nc;
          len_x = XDOM[xr].len;
          hi = len_x / nc_x;
          z = ZMAP[yr][xr];
          sa = ZON[z].st - ZON[z].ss;
          for (let i = 0; i < nc_x; i++) {
            abs_r[yr][xr] = abs_r[yr][xr] + sa * hj * hi * MFLUX[j_b][i_b];
            i_b = i_b + 1;
          }
        }
        j_b = j_b + 1;
      }
    }
    return abs_r;
  }

  // LEAKAGES
  else {
    let i_b, len_x, nc_x, hi, j_b, len_y, nc_y, hj;
    let M = QUAD['miu'].length;
    let leak = [0, 0, 0, 0];
    let ntc_x = 0,
      ntc_y = 0;
    for (let xr = 0; xr < XREG; xr++) {
      ntc_x = ntc_x + XDOM[xr].nc;
    }
    for (let yr = 0; yr < YREG; yr++) {
      ntc_y = ntc_y + YDOM[yr].nc;
    }
    j_b = 0;
    for (let yr = 0; yr < YREG; yr++) {
      nc_y = YDOM[yr].nc;
      len_y = YDOM[yr].len;
      hj = len_y / nc_y;
      for (let j = 0; j < nc_y; j++) {
        // LEFT
        for (let m = M / 4; m < (3 * M) / 4; m++) {
          let miu = QUAD['miu'][m],
            w = QUAD['w'][m];
          leak[0] = leak[0] + 0.5 * hj * w * miu * XFLOW[m][j_b][0];
        }
        // RIGHT
        for (let m = 0; m < M / 4; m++) {
          let miu = QUAD['miu'][m],
            w = QUAD['w'][m];
          leak[2] = leak[2] + 0.5 * hj * w * miu * XFLOW[m][j_b][ntc_x];
        }
        for (let m = (3 * M) / 4; m < M; m++) {
          let miu = QUAD['miu'][m],
            w = QUAD['w'][m];
          leak[2] = leak[2] + 0.5 * hj * w * miu * XFLOW[m][j_b][ntc_x];
        }
        j_b = j_b + 1;
      }
    }
    i_b = 0;
    for (let xr = 0; xr < XREG; xr++) {
      nc_x = XDOM[xr].nc;
      len_x = XDOM[xr].len;
      hi = len_x / nc_x;
      for (let i = 0; i < nc_x; i++) {
        // BOTTOM
        for (let m = M / 2; m < M; m++) {
          let theta = QUAD['theta'][m],
            w = QUAD['w'][m];
          leak[1] = leak[1] + 0.5 * hi * w * theta * YFLOW[m][0][i_b];
        }
        // TOP
        for (let m = 0; m < M / 2; m++) {
          let theta = QUAD['theta'][m],
            w = QUAD['w'][m];
          leak[3] = leak[3] + 0.5 * hi * w * theta * YFLOW[m][ntc_y][i_b];
        }
        i_b = i_b + 1;
      }
    }
    return leak;
  }
}

///////////////////////////////////////////////////////////////////////////// RESULTS SECTION => PLOT
/////////////////////////////////////////////////////////////////////////////////////////////////////
document.getElementById('results-section__plot--btn').onclick = function () {
  let options = document.querySelectorAll("input[name='plot_results']");
  options.forEach((option) => {
    if (option.checked) {
      selected = option.id;
      results_plot(selected);
    }
  });
};

// PLOT RESULTS
function results_plot(plot_option) {
  // VALIDATION
  if (results_validator() == false) {
    console_process_error('results');
    return;
  }

  document.getElementById('initial_img').style.display = 'none';
  if (document.getElementById('canvas_geom'))
    document.getElementById('canvas_geom').remove();
  if (document.getElementById('canvas_quad'))
    document.getElementById('canvas_quad').remove();
  if (document.getElementById('results__picture'))
    document.getElementById('results__picture').remove();

  let monitor = document.getElementById('main-content__layout--monitor');
  let plot = document.createElement('div');
  plot.id = 'results__picture';
  plot.style.width = '100%';
  plot.style.height = '100%';
  monitor.append(plot);
  let mydata = MFLUX;
  let plot_title = 'Scalar Flux';
  let z_axis_title = `\u03A6 [n cm${'-2'.sup()} s${'-1'.sup()}]`;
  if (plot_option == 'results-section__plot--absorption-rate') {
    mydata = results_data_absorption_rate();
    plot_title = 'Absorption Rate';
    z_axis_title = `abs [n s${'-1'.sup()}]`;
  }
  let plot_color = 'RdBu',
    rev = false;
  if (document.getElementById('results-section__plot--options__gray').checked) {
    plot_color = 'Greys';
    rev = false;
  }

  // 2D PLOT
  if (document.getElementById('results-section__plot--options__proy').checked) {
    Plotly.newPlot(
      'results__picture',
      [
        {
          // DATA
          x: results_x_layout(),
          y: results_y_layout(),
          z: mydata,
          type: 'contour',
          colorscale: plot_color,
          reversescale: rev,
          contours: {
            coloring: 'heatmap',
            showlabels: true,
            labelfont: {
              family: 'Raleway',
              size: 12,
              color: 'white',
            },
          },
          cmin: 0,
          colorbar: {
            title: z_axis_title,
            titleside: 'right',
            thickness: 10,
            len: 0.8,
          },
        },
      ],
      {
        // LAYOUT
        margin: {
          t: 50,
          b: 50,
        },
        xaxis: {title: 'X [cm]'},
        yaxis: {title: 'Y [cm]'},
        width: plot.offsetWidth,
        height: plot.offsetHeight,
      },
      {
        // CONFIG
        responsive: true,
        displaylogo: false,
        toImageButtonOptions: {
          format: 'png',
          filename: 'NTS_plot',
          height: 500,
          width: 700,
          scale: 1,
        },
      }
    );
  }

  // 3D PLOT
  else {
    Plotly.newPlot(
      'results__picture',
      [
        {
          // DATA
          x: results_x_layout(),
          y: results_y_layout(),
          z: mydata,
          type: 'surface',
          colorscale: plot_color,
          reversescale: rev,
          contours: {
            z: {
              show: true,
              usecolormap: true,
              highlightcolor: '#42f462',
              highlightwidth: 20,
              width: 20,
              start: 5,
              end: 10,
              project: {z: true},
            },
            coloring: 'heatmap',
            showlabels: true,
            labelfont: {
              family: 'Raleway',
              size: 12,
              color: 'white',
            },
          },
          cmin: 0,
          colorbar: {
            title: z_axis_title,
            titleside: 'right',
            thickness: 10,
            len: 0.8,
          },
        },
      ],
      {
        // LAYOUT
        margin: {
          t: 0,
          b: 0,
        },
        scene: {
          xaxis: {
            title: 'X [cm]',
            showgrid: true,
            showline: true,
            mirror: true,
            zeroline: false,
            gridcolor: 'rgb(240,240,240)',
          },
          yaxis: {
            title: 'Y [cm]',
            showgrid: true,
            showline: true,
            mirror: true,
            zeroline: false,
            gridcolor: 'rgb(240,240,240)',
          },
          zaxis: {
            title: z_axis_title,
            rangemode: 'tozero',
            autorange: true,
            showgrid: true,
            showline: true,
            mirror: true,
            showspikes: true,
            zeroline: false,
            gridcolor: 'rgb(240,240,240)',
          },
          camera: {
            type: 'orthographic',
            eye: {
              x: 1.75,
              y: 1.75,
              z: 1.75,
            },
          },
        },
        width: plot.offsetWidth,
        height: plot.offsetHeight,
      },
      {
        // CONFIG
        responsive: true,
        displaylogo: false,
        toImageButtonOptions: {
          format: 'png',
          filename: 'NTS_plot',
          height: 500,
          width: 700,
          scale: 1,
        },
      }
    );
  }
}

// X LAYOUT
function results_x_layout() {
  let xdata = [0];
  let x = 0;
  for (let xr = 0; xr < XREG; xr++) {
    let hi = XDOM[xr].len / XDOM[xr].nc;
    for (let i = 0; i < XDOM[xr].nc; i++) {
      x = x + hi;
      xdata.push(x);
    }
  }
  return xdata;
}

// Y LAYOUT
function results_y_layout() {
  let ydata = [0];
  let y = 0;
  for (let yr = 0; yr < YREG; yr++) {
    let hj = YDOM[yr].len / YDOM[yr].nc;
    for (let j = 0; j < YDOM[yr].nc; j++) {
      y = y + hj;
      ydata.push(y);
    }
  }
  return ydata;
}

// ABSORPTION RATE FOR PLOTTING
function results_data_absorption_rate() {
  let i_b, len_x, nc_x, hi, j_b, len_y, nc_y, hj, z, sa;
  let abs_data = JSON.parse(JSON.stringify(MFLUX));
  j_b = 0;
  for (let yr = 0; yr < YREG; yr++) {
    nc_y = YDOM[yr].nc;
    len_y = YDOM[yr].len;
    hj = len_y / nc_y;
    for (let j = 0; j < nc_y; j++) {
      i_b = 0;
      for (let xr = 0; xr < XREG; xr++) {
        nc_x = XDOM[xr].nc;
        len_x = XDOM[xr].len;
        hi = len_x / nc_x;
        z = ZMAP[yr][xr];
        sa = ZON[z].st - ZON[z].ss;
        for (let i = 0; i < nc_x; i++) {
          abs_data[j][i] = sa * hj * hi * MFLUX[j][i];
          i_b = i_b + 1;
        }
      }
      j_b = j_b + 1;
    }
  }
  return abs_data;
}

///////////////////////////////////////////////////////////////// RESULTS SECTION => DOWNLOAD RESULTS
/////////////////////////////////////////////////////////////////////////////////////////////////////
document.getElementById('results-section__download--btn').onclick =
  function () {
    if (results_validator()) {
      let storage = {
        ZN,
        ZON,
        XREG,
        XDOM,
        YREG,
        YDOM,
        ZMAP,
        QMAP,
        BC,
        QUAD,
        TOL,
        METH,
        ITER,
        CPU,
        MFLUX,
        MFLOW,
        XFLOW,
        YFLOW,
      };
      let dataBlob = new Blob([JSON.stringify(storage)], {
        type: 'application/json',
      });
      let url = URL.createObjectURL(dataBlob);
      this.setAttribute('href', url);
      this.setAttribute('download', 'NTS_results.json');
    } else {
      console_process_error('results');
    }
  };
