define(['app', 'app/model/library/author'], function (App) {

    return App.Collection.Authors = App.Backbone.Collection.extend({
        url: App.config.api_root + '/library/authors',
        model: App.Model.Author
    });

});
