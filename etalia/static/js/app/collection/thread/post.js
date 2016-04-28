define(['app', 'app/model/thread/post'], function (App) {

    return App.Collection.Posts = App.Backbone.Collection.extend({
        url: App.config.api_root + '/thread/posts',
        model: App.Model.Post
    });

});
