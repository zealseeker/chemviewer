const path = require('path');

module.exports = {
    entry: './src/main.js',
    output: {
        path: path.resolve(__dirname, 'chemviewer', 'static'),
    filename: 'bundle.js'
  },
  resolve: {
        alias: {
            'vue$': 'vue/dist/vue.esm.js'
        }
  },
  module: {
        rules:[
            {
                test: /\.(eot|svg|ttf|woff|woff2)(\?\S*)?$/,
                loader: 'file-loader'
            },
            {
                test: /\.css/i,
                use: ['style-loader', 'css-loader']
            }
        ]
  },
  mode: 'development',
  devtool: "eval-cheap-module-source-map"
};
