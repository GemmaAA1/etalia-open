define(['app', 'app/model/user/user'], function (App) {

    /*var defaults = {
        query: {
            view: 'nested'
        }
    };*/

    return App.Collection.Users = App.Backbone.PageableCollection.extend({
        url: App.config.api_root + '/user/users',
        //mode: 'infinite',
        model: App.Model.User,

        initialize: function() {
            //options = App.defaults(defaults, options);
            //this.queryParams = options.query;
        }
    });

});
