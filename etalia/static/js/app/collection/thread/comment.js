define(['app', 'app/model/thread/comment'], function (App) {

    return App.Collection.Comments = App.Backbone.Collection.extend({
        url: App.config.api_root + '/thread/comments',
        model: App.Model.Comment
    });

});
