let crtpath = '.'
const filemaxsize = 100  // K
ELEMENT.locale(ELEMENT.lang.en)
var app = new Vue({
  el: '#app',
  data: {
      text_smiles:'c1ccccc1 beneze\nCCC dimethylmethane',
      path_input: crtpath,
      filenames: [],
      dirs: [],
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
      inputShow: true
  },
  methods: {
    select_dir: function(dir){
      if (dir == 0){
        this.path_input = crtpath.replace(/\/[\u4e00-\u9fa5_a-zA-Z0-9]+$/,'')
      }else{
        this.path_input = crtpath + '/' + dir
      }
      this.list_dir()
    },
    list_dir: function(){
      // path = this.path_input.replace(/\./g,"%26")
      fetch('/dir?path='+this.path_input).then(res=>{
        return res.json()
      })
      .then(res=>{
        crtpath = this.path_input.replace(/\/$/,'')
        this.filenames=res['file']
        this.dirs = res['dir']
      })
    },
    select_file: function(filepath, abs){
      if(abs){
        path = filepath
      }else{
        path = crtpath+'/'+filepath
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
          this.loading=false
          return this.parseRes(res)
        }
      })
    },
    openDialog: function(cmp){
      console.log(cmp)
    },
    submit: function(){
      files = $('#upload_input')[0].files
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
        this.inputShow=false
        return this.parseRes(res)
      })
      .catch(err=>{
        this.loading=false
      })
    },
    parseRes(res){
      if (res && res.tableData){
        headers = []
        //TODO this can be improved
        this.isCompound = false
        this.isReaction = false
        if (res.type=='compound'){
          this.isCompound = true
        }else if(res.type=='reaction'){
          this.isReaction = true
        }
        for (i in res.tableData[0]){
          if (i!='_svg') headers.push(i)
        }
        this.allowed_cols = headers
        this.cols = headers
      }
      this.tableData = res.tableData.slice(0,100)
      this.allData = res.tableData
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
        template = this.get_svg_template(150, 100)
        for (i in this.tableData){
          item = this.tableData[i]
          if (item._svg){
            item._svg=item._svg.replace(/<svg[\s.\S]+?>/,template)
          }
        }
      }
    },
    get_svg_template(width, height, owidth, oheight){
      height = height || width
      owidth = owidth || 150
      oheight = oheight || owidth
      template = '<svg version="1.1" baseProfile="full" xmlns:svg="http://www.w3.org/2000/svg"' +
      'xmlns:rdkit="http://www.rdkit.org/xml" xmlns:xlink="http://www.w3.org/1999/xlink" ' +
      'xml:space="preserve" width="%width%px" height="%height%px" ' +
      'viewBox="0 0 %owidth% %oheight%">"'
      new_template = template.replace(/%width%/g,width).replace(/%height%/g,height)
        .replace(/%owidth%/g,owidth).replace(/%oheight%/,oheight)
      return new_template
    },
    change_pic_size(size){
      $('.grid_unit').css('width',size+'px')
      $('.grid_pic').css('width', size+'px')
      $('.grid_pic').css('height', size+'px')
      template = this.get_svg_template(size)
      for (i in this.tableData){
        item = this.tableData[i]
        if (item._svg){
          item._svg=item._svg.replace(/<svg[\s.\S]+?>/,template)
        }
      }
    },
    openDialog(cmp){
      this.dialogFormVisible=true
      if (this.isReaction){
        template = cmp._svg.match(/<svg[\s.\S]+?>/)[0]
        width = parseInt(template.match(/width='(\d*)px'/)[1])
        height = parseInt(template.match(/height='(\d*)px'/)[1])
        template = this.get_svg_template(900,300,width,height)
        $('.el-dialog').css('width', '920px')
      }else{
        template = this.get_svg_template(400,300)
        $('.el-dialog').css('width', '')
      }
      this.dialogPic = cmp._svg.replace(/<svg[\s.\S]+?>/, template)
      this.dialogCmp = cmp
    },
    handleCurrentChange(val) {
      console.log(val)
      this.tableData = this.allData.slice((val-1)*this.pageSize, val*this.pageSize)
    },
    click_input_hide(){
      this.inputShow = !this.inputShow
    }
  }
})
