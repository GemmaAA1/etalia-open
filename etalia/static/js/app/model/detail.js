define(['app'], function (App) {

    return App.Model.Detail = App.Backbone.Model.extend({
        defaults: {
            list: null,
            prev: null,
            next: null,
            view: null
        },

        initialize: function (options) {
            if (!this.get('list')) {
                throw '"list" is mandatory.';
            }
            if (!this.get('view')) {
                throw '"view" is mandatory.';
            }
        }
    });

});
