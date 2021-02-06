import Vue from 'vue';
import ElementUI from 'element-ui';
import 'element-ui/lib/theme-chalk/index.css';
import locale from 'element-ui/lib/locale/lang/en'
import {filemixin} from './filesys'
import {scatter_plot} from './d3fnc'

Vue.use(ElementUI, { locale })
Vue.config.devtools = true
const filemaxsize = 100  // K
function getUrlParam(name) {
  var reg = new RegExp("(^|&)" + name + "=([^&]*)(&|$)");
  var r = window.location.search.substr(1).match(reg);
  if (r != null) return decodeURI(r[2]); return null;
}
if (GLOBAL_PAGE == 'main'){
  var app = new Vue({
    el: '#app',
    mixins: [filemixin],
    data: {
        text_smiles:'c1ccccc1 beneze\nCCC dimethylmethane',
        activeTag: 'spreadsheet',
        tableData: [],
        allData: [],
        allowed_cols: ['SMILES', 'Name'],
        cols: ['SMILES', 'Name'],
        isCompound: false,
        isReaction: false,
        gridsize: 150,
        dialogCmp: null,
        dialogFormVisible: false,
        loading: false,
        pageSize: 100,
        inputShow: true,
        keepinput: false,
        filter_text: '',
        filters: [],
        filter_idx: [],
        dialogFilterVisible: false
    },
    methods: {
      select_file: function(filepath, abs){
        let path
        if(abs){
          path = filepath
        }else{
          path = this.crtpath+'/'+filepath
        }
        this.loading=true
        fetch('/file?path='+path)
        .then(res=>{
          if(res)
            return res.text()
        })
        .then(res=>{
          if(res){
            res = JSON.parse(res.replace(/\bNaN\b/g,null))
            if(!this.keepinput)
              this.inputShow=false
            this.loading=false
            return this.parseRes(res)
          }
        })
      },
      openDialog: function(cmp){
        console.log(cmp)
      },
      submit: function(){
        let files = $('#upload_input')[0].files
        var data = new FormData()
        if(files.length){
          // Upload file
          if(files[0].size > filemaxsize * 1024){
            alert('Too big file! Should be less than '+filemaxsize+'K')
            $('#upload_input')[0].value=''
            return false
          }
          data.append('file', files[0])
        }else{
          data.append('data', this.text_smiles)
        }
        this.loading=true
        fetch('/upload', {
          method: 'POST',
          body: data
        })
        .then(res=>{
          return res.text()
        })
        .then(res=>{
          res = JSON.parse(res.replace(/\bNaN\b/g,null))
          this.loading=false
          if(!this.keepinput)
            this.inputShow=false
          return this.parseRes(res)
        })
        .catch(err=>{
          this.loading=false
        })
      },
      parseRes(res){
        if (res && res.tableData){
          let headers = []
          //TODO this can be improved
          this.isCompound = false
          this.isReaction = false
          if (res.type=='compound'){
            this.isCompound = true
          }else if(res.type=='reaction'){
            this.isReaction = true
          }
          for (let i in res.tableData[0]){
            if (i!='_svg') headers.push(i)
          }
          this.allowed_cols = headers
          this.cols = headers
        }
        this.tableData = res.tableData.slice(0,100)
        this.allData = res.tableData
        this.filter_idx = Array.from({length: this.allData.length}, (_, d)=>d)
        if(res.tableData.length>100){
          $('#pagination').show()
        }else{
          $('#pagination').hide()
        }
        if(this.isReaction){
          this.handleTagClick({name: 'reaction'})
          this.activeTag = 'reaction'
        }else{
          this.handleTagClick({name: 'spreadsheet'})
          this.activeTag = 'spreadsheet'
        }
      },
      handleTagClick(tab, event){
        $('.tab_content').hide()
        $('#'+tab.name).show()
        if(tab.name=='grid'){
          this.gridsize=150
          this.change_pic_size(150)
        }else if(tab.name == 'spreadsheet'){
          let template = this.get_svg_template(150, 100)
          for (let i in this.tableData){
            let item = this.tableData[i]
            if (item._svg){
              item._svg=item._svg.replace(/<svg[\s.\S]+?>/,template)
            }
          }
        }
      },
      get_svg_template(width, height, owidth, oheight){
        height = height || width
        owidth = owidth || 200
        oheight = oheight || owidth
        let template = '<svg version="1.1" baseProfile="full" xmlns:svg="http://www.w3.org/2000/svg"' +
        'xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink" ' +
        'xml:space="preserve" width="%width%px" height="%height%px" ' +
        'viewBox="0 0 %owidth% %oheight%">"'
        let new_template = template.replace(/%width%/g,width).replace(/%height%/g,height)
          .replace(/%owidth%/g,owidth).replace(/%oheight%/,oheight)
        return new_template
      },
      change_pic_size(size){
        $('.grid_unit').css('width',size+'px')
        $('.grid_pic').css('width', size+'px')
        $('.grid_pic').css('height', size+'px')
        let template = this.get_svg_template(size)
        for (let i in this.tableData){
          let item = this.tableData[i]
          if (item._svg){
            item._svg=item._svg.replace(/<svg[\s.\S]+?>/,template)
          }
        }
      },
      openDialog(cmp){
        this.dialogFormVisible=true
        let template
        if (cmp._svg){
          if (this.isReaction){
            template = cmp._svg.match(/<svg[\s.\S]+?>/)[0]
            width = parseInt(template.match(/width='(\d*)px'/)[1])
            height = parseInt(template.match(/height='(\d*)px'/)[1])
            template = this.get_svg_template(900,300,width,height)
          }else{
            template = this.get_svg_template(400,300)
          }
          this.dialogPic = cmp._svg.replace(/<svg[\s.\S]+?>/, template)
        } else {
          this.dialogPic = 'Cannot show the picture'
        }
        this.dialogCmp = cmp
      },
      handleCurrentChange(val) {
        let idx = this.filter_idx.slice((val-1)*this.pageSize, val*this.pageSize)
        this.tableData = []
        for(let i=0; i<idx.length; i++){
          this.tableData.push(this.allData[idx[i]])
        }
        this.handleTagClick({name: this.activeTag})
      },
      handleSizeChange(val){
        this.handleCurrentChange(1)
      },
      click_input_hide(){
        this.inputShow = !this.inputShow
      },
      show_filter_dialog(){
        this.dialogFilterVisible = true
      },
      filter_command(command){
        if(command == 'show_filter_dialog')
          this.show_filter_dialog()
        else if (command == 'clear_filter'){
          this.filter_text = ''
          this.run_filter(true)
        }
      },
      dialog_filter_ok(){
        this.run_filter()
        this.dialogFilterVisible = false
      },
      run_filter(clear){
        let res = this.filter_text.split('=')
        let col, value
        if (res.length == 2){
          col = res[0]
          value = res[1]
        } else if (!clear){
          return false;
        }
        this.filter_idx = []
        for(let i=0; i<this.allData.length; i++){
          if(clear || this.allData[i][col]==value){
            this.filter_idx.push(i)
          }
        }
        if(this.filter_idx.length>this.pageSize){
          $('#pagination').show()
        }else{
          $('#pagination').hide()
        }
        this.handleCurrentChange(1)
      }
    }
  })
  if (GLOBAL_STD){
    app.list_dir()
    let filepath = getUrlParam('file')
    if(filepath){
      app.select_file(filepath, true)
    }
  }
}
else if(GLOBAL_PAGE=='analysis'){
  var app = new Vue({
    el: "#app",
    mixins: [filemixin],
    data:{
      input_file: '',
      table_columns: [],
      loading: false,
      table_length: 0,
      smiles_column: null,
      fingerprint_types: ['ECFP4'],
      fingerprint: 'ECFP4',
      progress: 1,
      jobid: null,
      failed: false,
      cycle_size: 3.5,
      svgs: null,
      pca_data: [
        {'PCA_0': 0.01, 'PCA_1':0.2}, 
        {'PCA_0': 0.12, 'PCA_1':0.5},
        {'PCA_0': 0.123, 'PCA_1':-0.7},
        {'PCA_0': 1.2, 'PCA_1':-0.2}]
    },
    methods:{
      select_file(filename){
        let path
        path = this.crtpath+'/'+filename
        this.input_file = path
        this.loading=true
        fetch('/file?mode=pre&path='+path)
        .then(res=>{
          if(res)
            return res.text()
        })
        .then(res=>{
          if(res){
            res = JSON.parse(res)
            this.loading=false
            this.table_columns = res.columns
            this.table_length = res.length
            if (res.jobid){
              this.jobid= res.jobid
            }
          }
        })
      },
      calc_fingerprint(){
        let fp = this.fingerprint
        fetch('/analysis/calcfp?fp='+fp+'&smi_col='+this.smiles_column+'&path='+this.input_file)
        .then(res=>{
          if(res)
            return res.text()
        })
        .then(res=>{
          this.jobid = res
          console.log('Done! jobid is '+ this.jobid)
          this.get_progress()
        })
      },
      calc_pca(){
        let fp = this.fingerprint
        let url = '/analysis/pca?' + new URLSearchParams({
          jobid: this.jobid,
          fp: this.fingerprint
        })
        fetch(url).then(res=>res.text())
        .then(res=>{
          let data = JSON.parse(res)
          console.log('Get Data, lenth: '+ data.length)
          console.log(data[0])
          this.pca_data = data
          scatter_plot(data)
        })
      },
      load_pca(){
        let url='/analysis/pca?' + new URLSearchParams({
          jobid: this.jobid,
          load: true
        })
        fetch(url).then(res=>res.text())
        .then(res=>{
          let data = JSON.parse(res)
          scatter_plot.run(data, {clicked: this.scatter_clicked})
        })
      },
      scatter_test(){
        scatter_plot.run(this.pca_data)
      },
      refresh_scatter(){
        scatter_plot.change_size(this.cycle_size)
      },
      get_progress(callback){
        let url = '/progress?jobid='+this.jobid
        fetch(url).then(res=>res.text())
        .then(res=>{
          res = parseInt(res)
          if (res > 0 && res <= 100)
            this.progress = res
          else
            this.progress = 0
          if (this.failed){
            this.failed == false
            return ;
          }
          if (this.progress != 100){
            setTimeout(() => {
              this.get_progress(callback)
            }, 500);
          } else {
            if (callback)
              callback()
          }
        })
      },
      get_structures(){
        fetch('/analysis/draw-structure?' + new URLSearchParams({
          smi_col: this.smiles_column,
          jobid: this.jobid
        }))
        .then(res=>{
          if(res)
            return res.text()
        })
        .then(res=>{
          this.jobid = res
          console.log('Done! jobid is '+ this.jobid)
          this.get_progress(this.load_structures)
        })
      },
      load_structures(){
        fetch('/analysis/fetch-structure?jobid='+this.jobid)
        .then(res=>res.text())
        .then(res=>{
          this.svgs = JSON.parse(res)
          console.log("Load SVGs, length: "+this.svgs.length)
        })
      },
      scatter_clicked(d){
        console.log(d)
        if(this.svgs && this.svgs[d.index]){
          console.log('Find svg')
          $('#structure').html(this.svgs[d.index])
        }
      }
    }
  })
}
if (GLOBAL_STD){
  app.list_dir()
  let filepath = getUrlParam('file')
  if(filepath){
    app.select_file(filepath, true)
  }
}
