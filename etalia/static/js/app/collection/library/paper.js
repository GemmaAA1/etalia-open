define(['app', 'app/model/library/paper'], function (App) {

    return App.Collection.Papers = App.Backbone.Collection.extend({
        url: App.config.api_root + '/library/papers',
        model: App.Model.Paper
    });

});
