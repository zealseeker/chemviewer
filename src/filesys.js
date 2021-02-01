export const filemixin = {
    data: {
        crtpath: '.',
        path_input: '.',
        filenames: [],
        dirs: [],
    },
    methods:{
        select_dir(dir){
            if (dir == 0){
              this.path_input = this.crtpath.replace(/\/[\u4e00-\u9fa5_a-zA-Z0-9\s\-]+$/,'')
            }else{
              this.path_input = this.crtpath + '/' + dir
            }
            this.list_dir()
        },
        list_dir(){
            fetch('/dir?path='+this.path_input).then(res=>{
                return res.json()
            })
            .then(res=>{
                this.crtpath = this.path_input.replace(/\/$/,'')
                this.filenames=res['file']
                this.dirs = res['dir']
            })
        }
    }
}

